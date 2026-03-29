#!/usr/bin/python3

import sys
from graphics import *
from pydispatch import dispatcher
from time import sleep
from Tkinter import *
import keyboard
from bigfloat import *


loc = ("/home/P/Python/Chemia/oranz ph12.xlsx")

RozmiarX = 1366
RozmiarY = 704

MoveX = RozmiarX/2
MoveY = RozmiarY/2


win = GraphWin("one day...", RozmiarX, RozmiarY)
win.setBackground('white')

N = 5

def funkcja(x):
	if x == 0:
		return 0
	return x

def clear():
    for item in win.items[:]:
        item.undraw()
    win.update()

def wyswietl():
	clear()
	poprzedniX = -MoveX
	poprzedniY = funkcja(poprzedniX)
	nowyX = 0
	nowyY = 0
	for i in range(1, 683):
		nowyX = -MoveX + i*2
		nowyY = funkcja(nowyX)
		p1 = Point(MoveX + poprzedniX, MoveY - poprzedniY);
		p2 = Point(MoveX + nowyX, MoveY - nowyY)
		line = Line(p1, p2)
		line.draw(win)
		#line = Line(p1, p2)
		#line.draw(win)
		poprzedniX = nowyX
		poprzedniY = nowyY

	win.update()

def main():
	global MoveX
	global MoveY
	wyswietl()
	# punkt (RozmiarX/2, RozmiarY/2) to (0, 0)
	stary = Point(MoveX, MoveY)
	for i in range(0, 10000000):
		nowy = win.getMouse()
		if(nowy.getX() != stary.getX() or nowy.getY() != stary.getY()):
			MoveX = nowy.getX()
			MoveY = nowy.getY()
			wyswietl()
		stary = nowy
		time.sleep(0.01)
		

	win.close()

main()


print "Zakonczono pomyslnie!"