import odswriter as ods

for i in range(0, 5):
	with ods.writer(open("odees.ods","wb")) as odsfile:
		odsfile.writerow(["String", "ABCDEF123456", "123456"])
		odsfile.writerow(["Float", i, 123, 123.123])