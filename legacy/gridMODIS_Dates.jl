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
import JSON
# Parallel computing
using Distributed, SharedArrays
# Profiler
using Profile

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
    	"--Dict"
            help = "JSON dictionary file to use"
            arg_type = String
            default = "/home/cfranken/code/gitHub/Gridding/gridding/modis_all.json"
        "--outFile", "-o"
            help = "output filename"
            arg_type = String
            default = "MODIS_map.nc"
		"--monthly"
	        help = "Use time-steps in terms of months (not days)"
	        action = :store_true
        "--latMin"
            help = "Lower latitude bound"
            arg_type = Float32
            default = -90.0f0
        "--latMax"
            help = "Upper latitude bound"
            arg_type = Float32
            default = 90.0f0
        "--lonMin"
            help = "Lower longitude bound"
            arg_type = Float32
            default = -180.0f0
        "--lonMax"
            help = "Upper longitude bound"
            arg_type = Float32
            default = 180.0f0
        "--dLat"
            help = "latitude resolution"
            arg_type = Float32
            default = 1.0f0
        "--dLon"
            help = "longitude resolution"
            arg_type = Float32
            default = 1.0f0
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
	cT = length(startDate:dDay:stopDate)


    # Just lazy (too cumbersome in code as often used variables here)
    latMax = ar["latMax"]
    latMin = ar["latMin"]
    lonMax = ar["lonMax"]
    lonMin = ar["lonMin"]
    dLat = ar["dLat"]
    dLon = ar["dLon"]
    eps = dLat/100

    # Define spatial grid:
    lat = collect(latMin+dLat/2.:dLat:latMax-dLat/2.0+eps)
    lon = collect(lonMin+dLon/2.:dLon:lonMax-dLon/2.0+eps)
	println("Output file dimension (time/lon/lat):")
    println(cT, "/", length(lon),"/", length(lat))
    # Create output file:
    dsOut = Dataset(ar["outFile"],"c")
    defDim(dsOut,"lon",length(lon))
	defDim(dsOut,"lat",length(lat))
	defDim(dsOut,"time", cT)
	dsLat = defVar(dsOut,"lat",Float32,("lat",), attrib = ["units" => "degrees_north","long_name" => "Latitude"])
	dsLon = defVar(dsOut,"lon",Float32,("lon",), attrib = ["units" => "degrees_east","long_name" => "Longitude"])
	dsTime= defVar(dsOut,"time",Float32,("time",),attrib = ["units" => "days since 1970-01-01","long_name" => "Time (UTC), start of interval"])
	dsLat[:]=lat
	dsLon[:]=lon


	# Define a global attribute
	dsOut.attrib["title"] = "Awesome gridded file"

    # Define gridded variables:
    n=zeros(Float32,(length(lat),length(lon)))
    SIF = zeros(Float32,(length(lat),length(lon)))
    # Parse JSON files as dictionary
	jsonDict = JSON.parsefile(ar["Dict"])
    #d2 = jsonDict["basic"]
    dGrid = jsonDict["grid"]
	# Get file naming pattern (needs YYYY MM and DD in there)
	fPattern = jsonDict["filePattern"]
	# Get main folder for files:
	folder = jsonDict["folder"]

	NCDict= Dict{String, NCDatasets.CFVariable}()
	println("Creating NC datasets in output:")
	for (key, value) in dGrid
		print(key," ")
		NCDict[key] = defVar(dsOut,key,Float32,("time","lon","lat"),deflatelevel=4, fillvalue=-999)
	end
	println(" ")
	#dSIF = defVar(dsOut,"sif",Float32,("lon","lat"),deflatelevel=4, fillvalue=-999)
	dN = defVar(dsOut,"n",Float32,("time","lon","lat"),deflatelevel=4, fillvalue=-999)
    # Define data array
    mat_data= zeros(Float32,(length(lon),length(lat),1+length(dGrid)))


	# Loop through time:
	# Time counter
	cT = 1
	for d in startDate:dDay:stopDate

		files = String[];
		for di in d:Dates.Day(1):d+dDay-Dates.Day(1)

			#********* This needs to be updated to use the Day of Year (and not MM and DD)!! *********#
			filePattern = reduce(replace,["YYYY" => lpad(Dates.year(di),4,"0"), "MM" => lpad(Dates.month(di),2,"0"),  "DD" => lpad(Dates.day(di),2,"0")], init=fPattern)
			println(filePattern)
			files = [files;glob(filePattern, folder)]
		end
		#println(files)

    	# Loop through all files
	    for a in files
	        # Read NC file
	        fin = Dataset(a)
	        # Check lat/lon first to see what data to read in

			#********* Here you need to read in lat/lon from MODIS (table or calculate on the fly) *********#
			#lat_in = fin[d2["lat"]].var[:]
	        #lon_in = fin[d2["lon"]].var[:]
			# Call the variables lat_in and lon_in and then best
	        lat_in_ = fin[d2["lat_bnd"]].var[:]
	        lon_in_ = fin[d2["lon_bnd"]].var[:]
			dim = size(lat_in_)

			# Get indices within the lat/lon boudning box:
			idx = findall((minLat[:,1].>latMin).&(maxLat[:,1].<latMax).&(minLon[:,1].>lonMin).&(maxLon[:,1].<lonMax))

			# Read data only for non-empty indices
	        if length(idx) > 0
				#print(size(lat_in_))
	            mat_in =  zeros(Float32,(length(lat_in_[:,1]),length(dGrid)+1))
				dim = size(mat_in)
	            # Read in all entries defined in JSON file:
				co = 1
				# This should just read in the datasets:
	            for (key, value) in dGrid
					#println(key, value)
	            	mat_in[:,co]=fin[value].var[:]
					co += 1
	            end
	            mat_in[:,end].=1
				# This computes the indices into which the respective lats and lons are falling into.
	            iLat_ = ((lat_in[idx].-latMin)/(latMax-latMin)*length(lat)).+1
	            iLon_ = ((lon_in[idx].-lonMin)/(lonMax-lonMin)*length(lon)).+1
				# Once you have done this, we can chat about the gridding itself (just a few lines of code here)

	            println("Read ", a, " ", length(idx))
	        else
	            println("Read ", a, " ", length(idx))
	        end
	        close(fin)
	    end
		# Filter all data, set averages
		dims = size(mat_data)
		println("Averaging final product...")

		NN = mat_data[:,:,end]
		dN[cT,:,:]=NN
		# Write out time:
		dsTime[cT]=d
		co = 1
		for (key, value) in dGrid
			da = round.(mat_data[:,:,co]./mat_data[:,:,end],sigdigits=5)
			da[NN.<1e-10].=-999
			NCDict[key][cT,:,:]=da
			co += 1
		end
		cT += 1
		fill!(mat_data,0.0)
	end
	close(dsOut)
end

main()
