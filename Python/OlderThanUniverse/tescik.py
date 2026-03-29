
import sys
from graphics import *
from time import sleep
from Tkinter import *
import keyboard
from pynput import mouse
from pynput import keyboard

RozmiarX = 1366
RozmiarY = 704

MoveX = RozmiarX/2
MoveY = RozmiarY/2

win = GraphWin("one day...", RozmiarX, RozmiarY)
win.setBackground('white')

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
	listener.start()

def on_move(x, y):
    print('Pointer moved to {0}'.format(
        (x, y)))


def on_click(x, y, button, pressed):
    print('{0} at {1}'.format(
        'Pressed' if pressed else 'Released',
        (x, y)))
    #if not pressed:
        # Stop listener
    #    return False

def on_scroll(x, y, dx, dy):
    print('Scrolled {0} at {1}'.format(
        'down' if dy < 0 else 'up',
        (x, y)))

# Collect events until released
with mouse.Listener(
        on_move=on_move,
        on_click=on_click,
        on_scroll=on_scroll) as listener:
    listener.join()

def on_press(key):
    try:
        print('alphanumeric key {0} pressed'.format(
            key.char))
    except AttributeError:
        print('special key {0} pressed'.format(
            key))

def on_release(key):
    print('{0} released'.format(
        key))
    if key == keyboard.Key.esc:
        # Stop listener
        return False

# Collect events until released
with keyboard.Listener(
        on_press=on_press,
        on_release=on_release) as listener:
    listener.join()

# ...or, in a non-blocking fashion:
#listener = mouse.Listener(
#    on_press=on_press,
#    on_release=on_release)

main()