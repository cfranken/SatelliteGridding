#!/home/cfranken//julia

# Argument Parser
using ArgParse
using Base, Dates, Printf
# NetCDF tools for reading and writing
using NCDatasets
# Basic statistics
using Statistics
# File search and completion
using Glob
# JSON files
using JSON
# Parallel computing
#using Distributed, SharedArrays
# Profiler
#using Profile

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--Dict"
            help = "JSON dictionary file to use"
            arg_type = String
            default = "/home/cfranken/code/gitHub/Gridding/gridding/tropomi_all.json"
        
        "--monthly"
            help = "Use time-steps in terms of months (not days)"
            action = :store_true
        "--startDate"
                help = "Start Date (in YYYY-MM-DD)"
                arg_type = String
                default = "2018-03-07"
        "--stopDate"
                help = "Stop Date (in YYYY-MM-DD)"
                arg_type = String
                default = "2018-10-31"
        "--dDays"
                help = "Time steps in days (or months if --monthly is set)"
                arg_type = Int64
                default = 8
    end
    return parse_args(s)
end



# Still need to make sure the corners are read in properly!
function getNC_var(fin, path)

    loc = split(path ,r"/")
    #println(loc)
    if length(loc)==1
        return fin[path].var
    elseif length(loc)>1
        gr = []
        for i in 1:length(loc)-1
            if i==1
                gr = fin.group[loc[i]]
            else
                gr = gr.group[loc[i]]
            end
        end
        return gr[loc[end]].var
    end
    
end

function getNC_attrib(fin, path, attri)
    loc = split(path ,r"/")
    #println(loc)
    if length(loc)==1
        return fin[path].attrib[attri]
    elseif length(loc)>1
        gr = []
        for i in 1:length(loc)-1
            if i==1
                gr = fin.group[loc[i]]
            else
                gr = gr.group[loc[i]]
            end
        end
        return gr[loc[end]].attrib[attri]
    end
end


function main()

    #addprocs()
    # Parse command line arguments
    ar = parse_commandline()

    # Find files to be processed
    startDate = DateTime(ar["startDate"])
    stopDate = DateTime(ar["stopDate"])
    if ar["monthly"]
        dDay = Dates.Month(ar["dDays"])
    else
        dDay = Dates.Day(ar["dDays"])
    end
    println(startDate, " ", stopDate)
    
    jsonDict = JSON.parsefile(ar["Dict"])
    d2       = jsonDict["basic"]
    # Get file naming pattern (needs YYYY MM and DD in there)
    fPattern = jsonDict["filePattern"]
    # Get main folder for files:
    folder   = jsonDict["folder"]

    # Loop through time:
    # Time counter
    cT = 1
    for d in startDate:dDay:stopDate
        files = String[];
        for di in d:Dates.Day(1):d+dDay-Dates.Day(1)
            #println("$(@sprintf("%04i-%02i-%02i", Dates.year(di),Dates.month(di),Dates.day(di)))")

            filePattern = reduce(replace,["YYYY" => lpad(Dates.year(di),4,"0"), "MM" => lpad(Dates.month(di),2,"0"),  "DD" => lpad(Dates.day(di),2,"0")], init=fPattern)
            #println(filePattern, " ", folder)
            files = [files;glob(filePattern, folder)]
        end
        fileSize = Int[];
        for f in files
            fileSize = [fileSize;stat(f).size]
        end
        #println(files)

        # Loop through all files
        for a in files[fileSize.>0]

            try 
                fin = Dataset(a)
                lat_in_ = getNC_var(fin, d2["lat_bnd"])
		lat_in_[1,1,1,1]
		lon_in_ = getNC_var(fin, d2["lon_bnd"])
                lon_in_[1,1,1,1]
                println("Passed, ", a)
                close(fin)
            catch e
                println(e)
                println("Failed, ", a)
                rm(a)
            end
            
        end
    end
end

main()
