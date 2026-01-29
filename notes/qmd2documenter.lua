-- Pandoc Lua filter to convert Quarto QMD to Documenter.jl markdown
-- Usage: pandoc hamiltonian.qmd --from markdown --to markdown --lua-filter=qmd2documenter.lua --wrap=none

-- Store document title from YAML frontmatter
local doc_title = nil

-- Map of section IDs to their titles (for @sec-name references)
local section_titles = {}

-- Macro replacements for Unicode output (adapted from expand_macros.lua)
-- Note: Pandoc expands \newcommand macros, so we match the expanded forms.

-- Hilbert space macros: \mathcal{H}_{...} or \mathcal{I}_{...} -> Unicode
-- The pattern includes the \! negative space that's in the original definitions
-- Note: Pandoc doubles the backslash in the internal representation (\\! instead of \!)
local hilbert_expanded = {
    -- These match the patterns from \newcommand definitions in the QMD header
    -- With doubled backslash as Pandoc stores them
    ["\\mathcal{H}_{\\\\!S}"] = "ℋ_S",
    ["\\mathcal{H}_{\\\\!I}"] = "ℋ_I",
    ["\\mathcal{H}_{\\\\!O}"] = "ℋ_O",
    ["\\mathcal{H}_{\\\\!GE}"] = "ℋ_{GE}",
    ["\\mathcal{H}_{\\\\!OS}"] = "ℋ_{OS}",
    ["\\mathcal{I}_{\\\\!M}"] = "ℐ_M",
    ["\\mathcal{I}_{\\\\!G}"] = "ℐ_G",
    ["\\mathcal{I}_{\\\\!E}"] = "ℐ_E",
    -- Also with single backslash (in case of different processing)
    ["\\mathcal{H}_{\\!S}"] = "ℋ_S",
    ["\\mathcal{H}_{\\!I}"] = "ℋ_I",
    ["\\mathcal{H}_{\\!O}"] = "ℋ_O",
    ["\\mathcal{H}_{\\!GE}"] = "ℋ_{GE}",
    ["\\mathcal{H}_{\\!OS}"] = "ℋ_{OS}",
    ["\\mathcal{I}_{\\!M}"] = "ℐ_M",
    ["\\mathcal{I}_{\\!G}"] = "ℐ_G",
    ["\\mathcal{I}_{\\!E}"] = "ℐ_E",
    -- Plain \mathcal{H} without subscript
    ["\\mathcal{H}"] = "ℋ",
}

-- Also handle un-expanded macros (in case Pandoc doesn't expand them)
local hilbert_macros = {
    ["Hilbert"] = "ℋ",
    ["HilbertS"] = "ℋ_S",
    ["HilbertI"] = "ℋ_I",
    ["HilbertO"] = "ℋ_O",
    ["HilbertGE"] = "ℋ_{GE}",
    ["HilbertOS"] = "ℋ_{OS}",
    ["HilbertM"] = "ℐ_M",
    ["HilbertG"] = "ℐ_G",
    ["HilbertE"] = "ℐ_E",
}

local operator_macros = {
    ["RWA"] = "\\text{RWA}",
    ["lab"] = "\\text{lab}",
    ["diag"] = "\\text{diag}",
    ["tr"] = "\\operatorname{tr}",
}

-- Also match the expanded forms from \textrm definitions
local operator_expanded = {
    ["\\textrm{tr}"] = "\\operatorname{tr}",
    ["\\textrm{diag}"] = "\\operatorname{diag}",
    ["\\textrm{RWA}"] = "\\text{RWA}",
    ["\\textrm{lab}"] = "\\text{lab}",
}

-- Find the matching closing brace for an opening brace at the given position
local function find_matching_brace(str, start_pos)
    local brace_count = 1
    local i = start_pos + 1
    while i <= #str do
        local char = str:sub(i, i)
        if char == '{' then
            brace_count = brace_count + 1
        elseif char == '}' then
            brace_count = brace_count - 1
            if brace_count == 0 then
                return i
            end
        end
        i = i + 1
    end
    return nil
end

-- Escape special characters for use in Lua pattern
local function escape_pattern(s)
    return s:gsub("([%^%$%(%)%%%.%[%]%*%+%-%?])", "%%%1")
end

-- Replace macros in math expressions and convert to Unicode where appropriate
local function expand_macros(math_str)
    -- Handle already-expanded Hilbert space macros -> Unicode
    -- Use plain string replacement with gsub (with escape_pattern)
    for pattern, replacement in pairs(hilbert_expanded) do
        -- Use plain text replacement by escaping pattern special chars
        local escaped = escape_pattern(pattern)
        math_str = math_str:gsub(escaped, replacement)
    end

    -- Remove any remaining \! (negative space) that might be left over
    math_str = math_str:gsub("\\\\!", "")
    math_str = math_str:gsub("\\!", "")

    -- Simplify subscripts: _{X} -> _X for single characters (matches target format)
    math_str = math_str:gsub("_{([%a])}", "_%1")

    -- Convert arrow macros to Unicode
    math_str = math_str:gsub("\\uparrow", "↑")
    math_str = math_str:gsub("\\downarrow", "↓")

    -- Handle un-expanded Hilbert space macros -> Unicode (in case Pandoc didn't expand them)
    for macro, replacement in pairs(hilbert_macros) do
        -- Match the macro at word boundaries
        math_str = math_str:gsub("\\(" .. macro .. ")([^%w])", replacement .. "%2")
        math_str = math_str:gsub("\\(" .. macro .. ")$", replacement)
    end

    -- Handle already-expanded operator macros
    for pattern, replacement in pairs(operator_expanded) do
        local escaped = escape_pattern(pattern)
        math_str = math_str:gsub(escaped, replacement)
    end

    -- Handle \mathop{\mathrm{...}} -> \text{...} or \operatorname{...}
    math_str = math_str:gsub("\\mathop{\\mathrm{diag}}", "\\text{diag}")
    math_str = math_str:gsub("\\mathop{\\mathrm{tr}}", "\\operatorname{tr}")

    -- Handle un-expanded operator macros -> \text{} or \operatorname{}
    for macro, replacement in pairs(operator_macros) do
        math_str = math_str:gsub("\\(" .. macro .. ")([^%w])", replacement .. "%2")
        math_str = math_str:gsub("\\(" .. macro .. ")$", replacement)
    end

    -- Handle \ket{X} -> |X⟩
    local pattern = "\\ket{"
    local start_pos = math_str:find(pattern, 1, true)
    while start_pos do
        local content_start = start_pos + #pattern - 1
        local closing_pos = find_matching_brace(math_str, content_start)
        if closing_pos then
            local content = math_str:sub(content_start + 1, closing_pos - 1)
            local prefix = math_str:sub(1, start_pos - 1)
            local suffix = math_str:sub(closing_pos + 1)
            math_str = prefix .. "|" .. content .. "⟩" .. suffix
        else
            break
        end
        start_pos = math_str:find(pattern, 1, true)
    end

    -- Handle \bra{X} -> ⟨X|
    pattern = "\\bra{"
    start_pos = math_str:find(pattern, 1, true)
    while start_pos do
        local content_start = start_pos + #pattern - 1
        local closing_pos = find_matching_brace(math_str, content_start)
        if closing_pos then
            local content = math_str:sub(content_start + 1, closing_pos - 1)
            local prefix = math_str:sub(1, start_pos - 1)
            local suffix = math_str:sub(closing_pos + 1)
            math_str = prefix .. "⟨" .. content .. "|" .. suffix
        else
            break
        end
        start_pos = math_str:find(pattern, 1, true)
    end

    -- Handle \ketbra{X}{Y} -> |X⟩\!\bra{Y}
    pattern = "\\ketbra{"
    start_pos = math_str:find(pattern, 1, true)
    while start_pos do
        local first_arg_start = start_pos + #pattern - 1
        local first_closing_pos = find_matching_brace(math_str, first_arg_start)
        if first_closing_pos then
            local first_arg = math_str:sub(first_arg_start + 1, first_closing_pos - 1)
            local second_arg_start = first_closing_pos + 1
            if math_str:sub(second_arg_start, second_arg_start) == "{" then
                local second_closing_pos = find_matching_brace(math_str, second_arg_start)
                if second_closing_pos then
                    local second_arg = math_str:sub(second_arg_start + 1, second_closing_pos - 1)
                    local prefix = math_str:sub(1, start_pos - 1)
                    local suffix = math_str:sub(second_closing_pos + 1)
                    -- Use \! for negative space between ⟩ and \bra
                    math_str = prefix .. "|" .. first_arg .. "⟩\\!\\bra{" .. second_arg .. "}" .. suffix
                else
                    break
                end
            else
                break
            end
        else
            break
        end
        start_pos = math_str:find(pattern, 1, true)
    end

    -- Add spacing around binary operators for better readability
    -- Space before \otimes, \oplus when preceded by a non-space character
    math_str = math_str:gsub("([^%s])\\otimes", "%1 \\otimes")
    math_str = math_str:gsub("([^%s])\\oplus", "%1 \\oplus")
    -- Space around = when adjacent to Unicode Hilbert space symbols
    math_str = math_str:gsub("(ℋ_[%w{}]+)=", "%1 =")
    math_str = math_str:gsub("(ℐ_[%w{}]+)=", "%1 =")

    return math_str
end

-- Handle inline math: $...$ -> ``...``
-- We need to output as RawInline to get double backticks
function Math(el)
    if el.mathtype == "InlineMath" then
        local expanded = expand_macros(el.text)
        -- Use RawInline with format "markdown" to preserve double backticks
        return pandoc.RawInline("markdown", "``" .. expanded .. "``")
    end
    -- Display math is handled at the block level in Para
    return el
end

-- Handle paragraphs: convert display math paragraphs to code blocks
function Para(el)
    -- Check if this paragraph starts with DisplayMath
    if #el.content >= 1 and el.content[1].t == "Math" and el.content[1].mathtype == "DisplayMath" then
        local math_el = el.content[1]
        local text = math_el.text
        local label = nil

        -- Check for equation label in the remaining content
        -- Pattern: Math, Space, Str "{#eq-...}"
        if #el.content >= 3 and el.content[2].t == "Space" and el.content[3].t == "Str" then
            local label_str = el.content[3].text
            local extracted_label = label_str:match("^{#(eq%-[%w%-]+)}$")
            if extracted_label then
                label = extracted_label
            end
        end

        -- Expand macros
        text = expand_macros(text)

        -- Build the math block content
        local content
        if label then
            -- Wrap in \begin{equation}\label{...}...\end{equation}
            content = "\\begin{equation}\\label{" .. label .. "}\n" .. text .. "\n\\end{equation}"
        else
            content = text
        end

        -- Return as a RawBlock to get ```math without space (Documenter.jl format)
        return pandoc.RawBlock("markdown", "```math\n" .. content .. "\n```")
    end
    return el
end

-- Handle YAML metadata: extract title
function Meta(meta)
    if meta.title then
        doc_title = pandoc.utils.stringify(meta.title)
    end
    -- Clear all metadata to prevent YAML output
    return {}
end

-- First pass: collect section titles for cross-references
local function collect_headers(doc)
    for _, block in ipairs(doc.blocks) do
        if block.t == "Header" then
            local id = block.identifier
            if id and id:match("^sec%-") then
                local title = pandoc.utils.stringify(block.content)
                section_titles[id] = title
            end
        end
    end
end

-- Handle headers: remove {#sec-...} anchors and filter out References section
function Header(el)
    -- Store section title for cross-references
    local id = el.identifier
    if id and id:match("^sec%-") then
        local title = pandoc.utils.stringify(el.content)
        section_titles[id] = title
    end

    -- Remove "References" header (handled separately in Documenter)
    local title = pandoc.utils.stringify(el.content)
    if title == "References" then
        return {}
    end

    -- Clear the identifier (Documenter doesn't use custom IDs)
    el.identifier = ""
    el.classes = {}
    el.attributes = {}

    return el
end

-- Handle citations and cross-references
-- Pandoc parses @sec-name, @eq-name, @fig-name, and [@key] all as Cite elements
function Cite(el)
    -- Check if this is a single citation (cross-reference)
    if #el.citations == 1 then
        local id = el.citations[1].id

        -- Handle @sec-name references -> [Section Title](@ref)
        if id:match("^sec%-") then
            local title = section_titles[id]
            if title then
                return pandoc.Link(pandoc.Str(title), "@ref")
            else
                -- Fallback: use the ID without sec- prefix
                return pandoc.Link(pandoc.Str(id), "@ref")
            end
        end

        -- Handle @eq-name references -> Equation ``\eqref{eq-name}``
        if id:match("^eq%-") then
            return pandoc.Inlines({
                pandoc.Str("Equation "),
                pandoc.RawInline("markdown", "``\\eqref{" .. id .. "}``")
            })
        end

        -- Handle @fig-name references -> "the figure above"
        if id:match("^fig%-") then
            return pandoc.Str("the figure above")
        end
    end

    -- Handle regular citations: [@key1; @key2] -> [key1, key2](@cite)
    local keys = {}
    for _, citation in ipairs(el.citations) do
        table.insert(keys, citation.id)
    end
    local cite_text = table.concat(keys, ", ")
    return pandoc.Link(pandoc.Str(cite_text), "@cite")
end

-- Handle images: convert paths from images/*.pdf to assets/*.svg
function Image(el)
    local src = el.src
    -- Convert images/X.pdf to assets/X.svg
    if src:match("^images/.*%.pdf$") then
        src = src:gsub("^images/", "assets/")
        src = src:gsub("%.pdf$", ".svg")
        el.src = src
    end
    return el
end

-- Handle Figure elements: remove identifier and convert nested Image paths
function Figure(el)
    -- Remove figure identifier (Documenter doesn't use {#fig-...})
    el.identifier = ""
    el.classes = {}
    el.attributes = {}

    -- Walk through content to convert image paths
    el.content = pandoc.walk_block(pandoc.Div(el.content), {
        Image = function(img)
            local src = img.src
            if src:match("^images/.*%.pdf$") then
                src = src:gsub("^images/", "assets/")
                src = src:gsub("%.pdf$", ".svg")
                img.src = src
            end
            return img
        end
    }).content

    return el
end

-- Handle emphasis: use underscores instead of asterisks
function Emph(el)
    -- Return raw markdown with underscores
    local content = pandoc.utils.stringify(el)
    return pandoc.RawInline("markdown", "_" .. content .. "_")
end

-- Handle divs (like ::: {#refs} :::) - remove them
function Div(el)
    -- Remove reference divs and other special Quarto divs
    if el.identifier == "refs" then
        return {}
    end
    return el
end

-- Handle inline sequences to fix edge cases like `Ω`$_{-}$ -> Ω₋
-- This walks through inline elements and merges patterns
function Inlines(inlines)
    local result = pandoc.Inlines({})
    local i = 1
    while i <= #inlines do
        local el = inlines[i]
        -- Check for pattern: Code followed by RawInline with subscript
        if el.t == "Code" and i + 1 <= #inlines then
            local next = inlines[i + 1]
            if next.t == "RawInline" and next.format == "markdown" then
                local subscript = next.text:match("^``_{([%-%+])}``$")
                if subscript then
                    -- Merge: Code + subscript -> Code with Unicode subscript
                    local unicodeSub = subscript == "-" and "₋" or "₊"
                    result:insert(pandoc.Code(el.text .. unicodeSub))
                    i = i + 2  -- Skip both elements
                else
                    result:insert(el)
                    i = i + 1
                end
            else
                result:insert(el)
                i = i + 1
            end
        else
            result:insert(el)
            i = i + 1
        end
    end
    return result
end

-- Main filter: process document
function Pandoc(doc)
    -- First pass: collect section titles
    collect_headers(doc)

    -- Process all blocks
    local new_blocks = {}

    -- Add title as H1 if we have one
    if doc_title then
        table.insert(new_blocks, pandoc.Header(1, pandoc.Str(doc_title)))
    end

    -- Process remaining blocks
    for _, block in ipairs(doc.blocks) do
        table.insert(new_blocks, block)
    end

    doc.blocks = new_blocks

    -- Return document with empty metadata
    return pandoc.Pandoc(doc.blocks, pandoc.Meta({}))
end

-- Return filters in order of execution
return {
    {Meta = Meta},
    {Pandoc = Pandoc},
    {Header = Header},
    {Para = Para},
    {Math = Math},
    {Cite = Cite},
    {Image = Image},
    {Figure = Figure},
    {Div = Div},
    {Emph = Emph},
    {Inlines = Inlines},
}
