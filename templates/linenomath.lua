-- quick fix to make line numbers work with `align` environment
function RawInline(e)
    if e.text:match "\\begin{align}" or e.text:match "\\begin{equation}" then
        return pandoc.RawInline(e.format, "\\begin{linenomath}"..e.text.."\\end{linenomath}")
    end
end
