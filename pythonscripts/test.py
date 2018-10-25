ids_obj = {}
with open("test") as input_file:
	input_file.readline()	
	for line in input_file:
		word = line.strip().split("\t")
		description = word[1]		
		try:		
			ids = word[4]
			ids = ids.split("/")
			for id in ids:
				if id in ids_obj:
					ids_obj[id].append(description)
				elif id not in ids_obj:
					ids_obj[id] = [description]
		except IndexError:
			pass
			
with open("promogeneinfo") as read_file:
	header = read_file.readline().strip()
	print(header+"\t"+"FUNCTIONS")
	for line in read_file:
		word = line.split("\t")
		entrez = word[0]
		try:
			functions = ",".join(ids_obj[entrez])
		except KeyError:
			functions = " "
		print(line.strip()+"\t"+functions)
		
