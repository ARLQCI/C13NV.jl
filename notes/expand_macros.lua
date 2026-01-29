-- List of macros (zero args) to replace, with replacement
local macro0_replacements = {
    {macroname = "tr", replace = "\\textrm{tr}"},
    {macroname = "diag", replace = "\\textrm{diag}"},
    {macroname = "RWA", replace = "\\textrm{RWA}"},
    {macroname = "lab", replace = "\\textrm{lab}"},
    {macroname = "Hilbert", replace = "\\mathcal{H}"},
    {macroname = "HilbertS", replace = "\\mathcal{H}_{\\!S}"},
    {macroname = "HilbertI", replace = "\\mathcal{H}_{\\!I}"},
    {macroname = "HilbertGE", replace = "\\mathcal{H}_{\\!GE}"},
    {macroname = "HilbertOS", replace = "\\mathcal{H}_{\\!OS}"},
    {macroname = "HilbertO", replace = "\\mathcal{H}_{\\!O}"},
    {macroname = "HilbertM", replace = "\\mathcal{I}_{\\!M}"},
    {macroname = "HilbertG", replace = "\\mathcal{I}_{\\!G}"},
    {macroname = "HilbertE", replace = "\\mathcal{I}_{\\!E}"},
}

-- List of macros (one arg) to replace, each with pre and post text
local macro1_replacements = {
    {macroname = "ket", replace_pre = "\\vert ", replace_post = " \\rangle"},
    {macroname = "bra", replace_pre = "\\langle ", replace_post = " \\vert"}
}

-- List of macros (two args) to replace, each with pre, mid, and post text
local macro2_replacements = {
    {macroname = "ketbra", replace_pre = "\\vert ", replace_mid = " \\rangle\\!\\langle ", replace_post = " \\vert"},
}

-- Find the matching closing brace for an opening brace at the given position
function find_matching_brace(str, start_pos)
    local brace_count = 1
    local i = start_pos + 1
    while i <= #str do
        local char = str:sub(i,i)
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
    return nil  -- no matching brace found
end

-- Replace macros in math expressions
function replace_macros_in_math(math_str)

    -- Handle all macros from the zero-argument replacement list
    for _, macro in ipairs(macro0_replacements) do
        -- patterns match if followed by the non-word character
        -- we use a "frontier pattern" (%f) so that it also matches at "end of word"
        local pattern = "\\" .. macro.macroname .. "%f[%W]"
        math_str = math_str:gsub(pattern, macro.replace)
    end

    -- Handle all macros from the single argument replacement list
    for _, macro in ipairs(macro1_replacements) do
        local pattern = "\\" .. macro.macroname .. "{"
        local start_pos = math_str:find(pattern)

        -- Keep replacing while we find instances of the macro
        while start_pos do
            local content_start = start_pos + #pattern - 1
            local closing_pos = find_matching_brace(math_str, content_start)

            if closing_pos then
                -- Extract content between braces
                local content = math_str:sub(content_start + 1, closing_pos - 1)

                -- Replace macro with pre + content + post
                local prefix = math_str:sub(1, start_pos - 1)
                local suffix = math_str:sub(closing_pos + 1)
                math_str = prefix .. macro.replace_pre .. content .. macro.replace_post .. suffix
            else
                -- No matching brace found - skip this occurrence
                break
            end

            -- Look for next occurrence
            start_pos = math_str:find(pattern)
        end
    end

    -- Handle macros from the two-argument replacement list
    for _, macro in ipairs(macro2_replacements) do
        local pattern = "\\" .. macro.macroname .. "{"
        local start_pos = math_str:find(pattern)

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
                        math_str = prefix .. macro.replace_pre .. first_arg .. macro.replace_mid .. second_arg .. macro.replace_post .. suffix
                    else
                        break
                    end
                else
                    break
                end
            else
                break
            end

            start_pos = math_str:find(pattern)
        end
    end

    return math_str

end


function Math(el)
    el.text = replace_macros_in_math(el.text)
    return el
end


function DisplayMath(el)
    el.text = replace_macros_in_math(el.text)
    return el
end
