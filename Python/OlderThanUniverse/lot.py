#!/usr/bin/python3

import sys
from graphics import *
from time import sleep
from Tkinter import *
import keyboard
from math import *

RozmiarX = 1366
RozmiarY = 704

sx = 0
sy = 0

rot = 0

vx = 0
vy = 5

dx = 200
dy = 0

A = 1
B = 1
W = 5

MoveX = RozmiarX/2
MoveY = RozmiarY/2

win = GraphWin("trajektoria", RozmiarX, RozmiarY)
win.setBackground('white')

def wyswietl(numer):
	global sx, sy, vx, vy, rot, dx, dy, A, B, W
	# tutaj jest obliczanie
	drot = atan2(dy - sy, dx - sx)
	vrot = atan2(vy, vx)

	V2 = vx*vx + vy*vy
	D = sqrt((dx - sx)*(dx - sx) + (dy - sy)*(dy - sy))
	C = erf((1 - cos(drot - vrot)*V2/(2*A*D))/W)

	Crot = (drot - vrot)/2 + C*(drot + vrot)/2

	if Crot > pi:
		Crot -= 2*pi
	if Crot < -pi:
		Crot += 2*pi
	#print(rot, "   ", drot)
	if rot < Crot:
		if rot < Crot - pi:
			rot -= (rot + 2*pi - Crot)/B
		else:
			rot += (Crot - rot)/B
	else:
		if rot > Crot + pi:
			rot += (Crot + 2*pi - rot)/B
		else:
			rot -= (rot - Crot)/B

	if rot > pi:
		rot -= 2*pi
	if rot < -pi:
		rot += 2*pi

	print(C, '    ', rot, '   ', Crot)
	vx += A*C * cos(rot)
	vy += A*C * sin(rot)

	#if V2 > A*D*2:
	#	print("---")
	#	vx -= A * cos(rot) 
	#	vy -= A * sin(rot)
	#else:
	#	print("+++")
	#	vx += A * cos(rot)
	#	vy += A * sin(rot)
	
	sx += vx
	sy += vy
	# koniec
	p = Point(MoveX + sx, MoveY - sy)
	p.draw(win)

	dp = Point(MoveX + dx, MoveY - dy)
	dp.draw(win)

	win.update()

def main():
	global MoveX
	global MoveY

	for i in range(0, 10000000):
		if i > 10:
			wyswietl(i - 10)
		time.sleep(0.1)

	win.close()

main()


print "Zakonczono pomyslnie!"