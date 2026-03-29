#!/usr/bin/python3

import sys

plik = open('Chemia/test.txt')
s = plik.read()

tab1 = [0 for i in range(0, 3)]
tab2 = [0 for i in range(0, 3)]
tab3 = [0 for i in range(0, 3)]
tab4 = [0 for i in range(0, 3)]

def convert():
	global tab1, tab2, tab3, tab4, s
	number = 1
	ind = 0
	dot = -1
	x = 0
	START = False
	for i in range(0, len(s)):
		if 0 <= ord(s[i]) - 48 and ord(s[i]) - 48 <= 9:
			if dot >= 0:
				dot = dot + 1
			if START:
				x *= 10
				x += ord(s[i]) - 48
			else:
				START = True
				x = ord(s[i]) - 48
		elif s[i] == ',':
			dot = 0
		elif START:
			START = False
			if dot >= 0:
				x *= 1.0000000000001
				x /= 10**dot
			if number == 1:
				tab1[ind] = x
			elif number == 2:
				tab2[ind] = x
			elif number == 3:
				tab3[ind] = x
			else:
				tab4[ind] = x
				number = 0
				ind = ind + 1
			number = number + 1
			x = 0
			dot = -1

convert()

for i in range(0, 3):
	print tab3[i]
	
print "Zakonczono pomyslnie!"