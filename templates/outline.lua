str = pandoc.utils.stringify

-- remove headers with class `outline` in latex output
function Header(e)
    if e.classes:includes("outline") and FORMAT:match "latex" then
        return {}
    end
end

-- Check if paragraph contains separator `––` either at the beginning (outline
-- without content) or in the middle (outline with content). For LaTeX output,
-- the outline is removed, for others it is placed in a separate element.
function Para(e)

    local elems = e.content
    local m, i = elems:find_if(is_marker)
    if m == nil then
        return e
    end

    --print("current:", str(elems):sub(1,20))

    -- handle empty outlines
    if i == 1 and elems[2].tag == "Space" then
        --print("found an empty outline", #elems)
        elems:remove(1)
        elems:remove(1)
        if FORMAT:match "latex" then
            return {}
        else
            return pandoc.Div({ parse_outline(elems) }, {class="outlined empty"})
        end
    end

    -- handle outline with paragraph
    if elems[i-1].tag == "SoftBreak" and elems[i+1].tag == "SoftBreak" then
        --print("found an outline with paragraph", #elems)
        local outline = {}
        for k = 1,i-2 do
            outline[k] = elems:remove(1)
            --print("removing elem "..k, str(outline[k]))
        end
        elems:remove(1) -- SoftBreak
        elems:remove(1) -- divider
        elems:remove(1) -- SoftBreak
        outline = pandoc.List(outline)
        if FORMAT:match "latex" then
            return pandoc.Para(elems)
        else
            return pandoc.Div({ parse_outline(pandoc.List(outline)), pandoc.Para(elems) }, {class="outlined complete"})
        end
    end

    -- false positive
    return e
end

function is_marker(s)
    return s.tag == "Str" and s.text == "––"
end

function parse_outline(elems)
    local parts = {}
    local part_elems = {}
    for i,e in ipairs(elems) do
        if e.tag == "SoftBreak" and elems[i+1].tag == "Str" and elems[i+1].text == "-" and elems[i+2].tag == "Space" then
            parts[#parts+1] = pandoc.Plain(part_elems)
            part_elems = {}
        elseif #part_elems == 0 and e.tag == "Str" and e.text == "-" then
            -- do not add leading "-"
        elseif #part_elems == 0 and e.tag == "Space" then
            -- do not add leading space
        else
            part_elems[#part_elems+1] = e
        end
    end
    parts[#parts+1] = pandoc.Plain(part_elems)

    local outline = { pandoc.Div({ table.remove(parts, 1) }, { class="outline-primary" })  }
    if #parts > 0 then
        outline[2] = pandoc.BulletList(parts)
    end

    outline = pandoc.Div(outline, { class="outline" })
    return outline
end
