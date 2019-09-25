
#SampleID,ACAD Number,Sex Assignment,Location ,Region,Latitude,Longitude,Altitude,Site type,Bone Type,Collected From,Date,Endogenous,Endogenous2,Clonal,Raw_Mit,Uniq_Mit,Raw_Nuc,Uniq_Nuc,RatioMtNu

BEGIN {
	FS=","
	OFS="\t"
	print "sample", "site_type", "material", "material2", "collection", "lat", "lon", "alt", "age", "endog"
}

NR>1 {
	aid="A"$2
	#sex=$3
	site_type=$9
	material=$10
	collection=$11
	lat=$6
	lon=$7
	alt=$8
	age=$12
	endog=$13

	# Negative values are in Beringia and North America,
	# which are closer to Russia than to England, so make
	# the numbers reflect this.
	if (lon < -15)
		lon += 360

	if (site_type ~ /Mass death/)
		site_type = "Cave"
	if (site_type ~ /Anthropogenic/)
		# skip these
		next

	if (material ~ /Flat Bone/)
		material = "FlatBone"
	else if (material ~ /Long Bone/)
		material = "LongBone"
	else if (material ~ /Undefined/)
		material = "None"

	if (material ~ /Tooth|Crania/)
		material2 = "Crania"
	else if (material != "None")
		material2 = "Postcrania"
	else
		material2 = "None"


	print aid, site_type, material, material2, collection, lat, lon, alt, age, endog
}
