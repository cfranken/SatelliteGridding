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
using Distributed, SharedArrays
# Profiler
using Profile

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
    	"--Dict"
            help = "JSON dictionary file to use"
            arg_type = String
            default = "/home/cfranken/code/gitHub/Gridding/gridding/tropomi_all.json"
        "--outFile", "-o"
            help = "output filename (default OCO2_SIF_map.nc)"
            arg_type = String
            default = "OCO2_SIF_map.nc"
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

function getFilter(name, jsonDict)
	ff = try
		jsonDict[name]
	catch
		Dict()
	end
	return ff
end
# This splits up the entire region into one grid set
function divLine!(lat1,lon1,lat2,lon2,n, points,j )
    dLat = (lat2-lat1)/(2*n)
    dLon = (lon2-lon1)/(2*n)
	startLat = lat1+dLat
	startLon = lon1+dLon
	for i in 1:n
		#println(startLat+2*(i-1)*dLat, " ", startLon+2*(i-1)*dLon)
		#weights[(iLat-minLat+1), (iLon-minLon+1)]+=1
		points[j,i,1] = startLat+2*(i-1)*dLat
		points[j,i,2] = startLon+2*(i-1)*dLon
	end
end

# For first run (2 baselines)
function divLine2!(lat1,lon1,lat2,lon2,n, lats, lons)
    dLat = (lat2-lat1)/(2*n)
    dLon = (lon2-lon1)/(2*n)
	startLat = lat1+dLat
	startLon = lon1+dLon
	@inbounds for i in 1:n
	    lats[i] = startLat+2*(i-1)*dLat
	    lons[i] = startLon+2*(i-1)*dLon
	end
end

# Divide each polygon into multiple points
function getPoints!(points, vert_lat, vert_lon, n,lats_0, lons_0,lats_1, lons_1 )
	# Get reference points for two lines at the extremes:
	# println(vert_lat)
	divLine2!(vert_lat[1],vert_lon[1],vert_lat[2],vert_lon[2],n,lats_0, lons_0)
    divLine2!(vert_lat[4],vert_lon[4],vert_lat[3],vert_lon[3],n,lats_1, lons_1)
    for i in 1:n
         divLine!(lats_0[i], lons_0[i] ,lats_1[i], lons_1[i],n, points,i )
    end
end

function getNC_var(fin, path, DD::Bool)
	loc = split(path ,r"/")
	#println(loc)
	if length(loc)==1
		return fin[path].var[:]
	elseif length(loc)>1
		gr = []
		for i in 1:length(loc)-1
			if i==1
				gr = fin.group[loc[i]]
			else
				gr = gr.group[loc[i]]
			end
		end
		#println(loc[end])
		si = size(gr[loc[end]])
		# DD means there is a 2nd index for footprint bounds of dimension 4!
		if DD
			if si[1]==4
				return reshape(gr[loc[end]].var[:],4,prod(si[2:end]))'
			elseif si[end]==4
				return reshape(gr[loc[end]].var[:],prod(si[1:end-1]),4)
			end
		end
		return reshape(gr[loc[end]].var[:],prod(si))
		#a = reshape(gr[loc[end]],4,215*2906*1)
		#return gr[loc[end]].var[:]
	else
		println("Something is off in getNC_var")
	end

end

# Test for how fast we can compute point in distributed mode (not yet used)
function compPoints!(lat,lon,n)
	dim = size(lat)
	lats_0 = zeros(n)
	lons_0 = zeros(n)
	lats_1 = zeros(n)
	lons_1 = zeros(n)
	a = SharedArray{Float32}((dim[1],n,n,2))
	@distributed for i in 1:dim[1]
        po = getPoints!(a[i,:,:,:],lat[i,:],lon[i,:],n,lats_0, lons_0,lats_1, lons_1 )
    end
    return a
end

function favg_all!(arr,lat,lon,inp,s,s2,n, latMin, latMax, lonMin, lonMax, nLat, nLon, points)
	dim = size(arr)
	#println(dim)
	# Predefine some arrays to reduce allocations
	ix = zeros(Int32,n^2)
	iy = zeros(Int32,n^2)
	lats_0 = zeros(n)
	lons_0 = zeros(n)
	lats_1 = zeros(n)
	lons_1 = zeros(n)
	iLon = floor.(Int32,lon)
	iLat = floor.(Int32,lat)
	minLat = minimum(Int32,floor.(lat), dims=2)
	maxLat = maximum(Int32,floor.(lat), dims=2)
	minLon = minimum(Int32,floor.(lon), dims=2)
	maxLon = maximum(Int32,floor.(lon), dims=2)
	distLon = maxLon-minLon
	# How many individual grid cells might actually be there:
	global dimLat = maxLat-minLat
	global dimLon = maxLon-minLon
	fac = n^2
    @inbounds for i in 1:s
		#println(i, " ", dimLat[i], " ", dimLon[i])
		# Take it easy if all corners already fall into one grid box:
		if (dimLat[i]==1) & (dimLon[i]==1)
			for z in 1:s2
	                arr[iLon[i,1],iLat[i,1],z] +=  fac*inp[i,z]
	        end
		# if not, compute appropriate weights
		elseif (distLon[i])<n
	        getPoints!(points,lat[i,:],lon[i,:],n,lats_0, lons_0,lats_1, lons_1 )
	        ix[:] = floor.(Int32,points[:,:,1][:])
	        iy[:] = floor.(Int32,points[:,:,2][:])
			#iy[iy .> dim[1]].=dim[1]
			#iy[iy .< 1].=1

	        @inbounds for j in eachindex(ix)
	            for z in 1:s2
	                arr[iy[j],ix[j],z] +=  inp[i,z]
	                #println(arr[ix[j],iy[j],z])
	            end
	        end
		#else
			#println("Warning ", iLat[i], " ",iLon[i], " ", maxLon[i]-minLon[i])
		end
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
    d2 = jsonDict["basic"]
    dGrid = jsonDict["grid"]

	f_eq = getFilter("filter_eq",jsonDict)
	f_gt = getFilter("filter_gt",jsonDict)
	f_lt = getFilter("filter_lt",jsonDict)

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

	# Still hard-coded here, can be changed:
    nGrid = 10

    points = zeros(Float32,(nGrid,nGrid,2))
	#global indices = zeros(Int32,(nGrid,nGrid,2))
	global weights = zeros(Int32,(nGrid,nGrid))

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
	        # Read NC file
			#try
			println(a)
			#println(fStats.size)
	        fin = Dataset(a)
	        # Check lat/lon first to see what data to read in
	        #lat_in = fin[d2["lat"]].var[:]
	        #lon_in = fin[d2["lon"]].var[:]
            lat_a = getNC_var(fin, "GEOLOCATION/LatitudeA",true)
			lat_b = getNC_var(fin, "GEOLOCATION/LatitudeB",true)
			lat_c = getNC_var(fin, "GEOLOCATION/LatitudeC",true)
			lat_d = getNC_var(fin, "GEOLOCATION/LatitudeD",true)
			lon_a = getNC_var(fin, "GEOLOCATION/LongitudeA",true)
			lon_b = getNC_var(fin, "GEOLOCATION/LongitudeB",true)
			lon_c = getNC_var(fin, "GEOLOCATION/LongitudeC",true)
			lon_d = getNC_var(fin, "GEOLOCATION/LongitudeD",true)
			println(size(lon_d))
	        #lat_in_ = getNC_var(fin, d2["lat_bnd"],true)
	        #lon_in_ = getNC_var(fin, d2["lon_bnd"],true)
			lat_in_ = [lat_a lat_b lat_c lat_d]
	        lon_in_ = [lon_a lon_b lon_c lon_d]
			for i=1:length(lon_b)
				for j=1:4
				if lon_in_[i,j]>180
					lon_in_[i,j] = lon_in_[i,j]-360
				end
			end
			end
			println(size(lon_in_))
			dim = size(lat_in_)
			# Transpose if orders are swapped
			if dim[1]==4
				lat_in_ = lat_in_'
				lon_in_ = lon_in_'
			end
	        # Find all indices within lat/lon bounds:
			minLat = minimum(lat_in_, dims=2)
			maxLat = maximum(lat_in_, dims=2)
			minLon = minimum(lon_in_, dims=2)
			maxLon = maximum(lon_in_, dims=2)

			# Get indices within the lat/lon bounding box and check filter criteria:
			bool_add = (minLat[:,1].>latMin) .+ (maxLat[:,1].<latMax) .+ (minLon[:,1].>lonMin) .+ (maxLon[:,1].<lonMax) .+ ((maxLon[:,1].-minLon[:,1]).<50)
			bCounter = 5
			for (key, value) in f_eq

				bool_add += (getNC_var(fin, key,false).==value)
				bCounter+=1
			end
			for (key, value) in f_gt
				bool_add += (getNC_var(fin, key,false).>value)
				bCounter+=1
			end
			for (key, value) in f_lt
				bool_add += (getNC_var(fin, key,false).<value)
				bCounter+=1
			end
			#idx = findall((minLat[:,1].>latMin).&(maxLat[:,1].<latMax).&(minLon[:,1].>lonMin).&(maxLon[:,1].<lonMax))
			idx = findall(bool_add.==bCounter)
			#println(length(idx), " ", length(idx2))
			#println(bCounter)
			#println(fin[key].var[:])
			# Read data only for non-empty indices
	        if length(idx) > 0
				#print(size(lat_in_))
	            mat_in =  zeros(Float32,(length(lat_in_[:,1]),length(dGrid)+1))
				dim = size(mat_in)
	            # Read in all entries defined in JSON file:
				co = 1
	            for (key, value) in dGrid
					#println(key, value)
					mat_in[:,co]=getNC_var(fin, value,false)
					co += 1
	            end
	            mat_in[:,end].=1
	            iLat_ = ((lat_in_[idx,:].-latMin)/(latMax-latMin)*length(lat)).+1
	            iLon_ = ((lon_in_[idx,:].-lonMin)/(lonMax-lonMin)*length(lon)).+1
				@time  favg_all!(mat_data, iLat_,iLon_,mat_in[idx,:],length(idx),dim[2],nGrid, latMin, latMax, lonMin,lonMax, length(lat), length(lon), points )
	            println("Read ", a, " ", length(idx))
	        else
	            println("Read ", a, " ", length(idx))
	        end
	        close(fin)
			#catch
			#	println("Error in file caught")
			#end
	    end
		# Filter all data, set averages
		dims = size(mat_data)
		println("Averaging final product...")
		if maximum(mat_data[:,:,end])>0
			NN = mat_data[:,:,end]./nGrid^2
			dN[cT,:,:]=NN
			# Write out time:
			dsTime[cT]=d
			co = 1
			for (key, value) in dGrid
				da = round.(mat_data[:,:,co]./mat_data[:,:,end],sigdigits=5)
				da[NN.<1e-10].=-999
				@inbounds NCDict[key][cT,:,:]=da
				co += 1
			end
		else
			dN[cT,:,:]=0
			dsTime[cT]=d

		end
		cT += 1
		fill!(mat_data,0.0)
	end
	close(dsOut)
end

main()
