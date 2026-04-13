"""
    find_files(pattern::String, folder::String, date::DateTime) -> Vector{String}

Find input files matching the given pattern for a specific date.

Resolves `YYYY`, `YY`, `MM`, `DD`, and `DOY` placeholders in both `pattern` and `folder`,
then globs for matching files. Returns only non-empty files.
"""
function find_files(pattern::String, folder::String, date::DateTime)::Vector{String}
    file_pattern = _resolve_date_placeholders(pattern, date)
    folder_resolved = _resolve_date_placeholders(folder, date)
    files = glob(file_pattern, folder_resolved)
    filter(f -> stat(f).size > 0, files)
end

"""
    find_files_range(pattern, folder, start_date, end_date) -> Vector{String}

Find files for a date range (inclusive), one day at a time.
"""
function find_files_range(pattern::String, folder::String,
                          start_date::DateTime, end_date::DateTime)::Vector{String}
    files = String[]
    for d in start_date:Dates.Day(1):end_date
        append!(files, find_files(pattern, folder, d))
    end
    files
end

function _resolve_date_placeholders(s::String, date::DateTime)::String
    reduce(replace, [
        "YYYY" => lpad(Dates.year(date), 4, "0"),
        "YY"   => lpad(Dates.year(date) % 100, 2, "0"),
        "MM"   => lpad(Dates.month(date), 2, "0"),
        "DD"   => lpad(Dates.day(date), 2, "0"),
        "DOY"  => lpad(Dates.dayofyear(date), 3, "0"),
    ]; init=s)
end
