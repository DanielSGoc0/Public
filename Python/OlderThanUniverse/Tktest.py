#!/usr/bin/python3

from Tkinter import *
import Tkinter


x = 0

root = Tk()
frame = Frame(root)
frame.pack()

cord = [0, 0, 0, 0]

C = Tkinter.Canvas(root, height=670, width=1366, bg = 'white')
C.pack()

line_old = C.create_line(cord, fill = 'white')
line_new = C.create_line(cord)

def funkcja(i):
	global x
	if i % 2 == 0:
		return 10*i
	else:
		return (i - x) * (i - x)/10

def Action1():
	global x
	repaint()
	x = x + 1
	print x
	paint()

def Action2():
	global x
	repaint()
	x = x - 1
	print x
	paint()

def repaint():
	global x, line_old, cord
	line_old = C.create_line(cord, fill = 'white')

def paint():
	global x, line_new, cord
	cord = [funkcja(i) for i in range(0, 100)]
	line_new = C.create_line(cord)


Button1 = Tkinter.Button(frame, text ="to the right", command = Action1)
Button1.pack(side = RIGHT)

Button2 = Tkinter.Button(frame, text ="to the left", command = Action2)
Button2.pack(side = RIGHT)

root.mainloop()