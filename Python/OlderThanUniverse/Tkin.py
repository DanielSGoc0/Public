#!/usr/bin/python3

from Tkinter import *

RozmiarX = 1366
RozmiarY = 704

class Window(Frame):

	def __init__(self, master=None):
		Frame.__init__(self, master)
		self.master = master
		self.create()

	def create(self):
		self.master.title("python tkinter go!")
		self.pack(fill=BOTH, expand=1)
		quitButton = Button(self, text="exit!!!", command=self.exit_program)
		quitButton.place(x=RozmiarX/2, y=RozmiarY/2)

		pasek = Menu(self.master)
		self.master.config(menu=pasek)

		file = Menu(pasek)
		file.add_command(label="exit!", command=self.exit_program)
		pasek.add_cascade(label="plik", menu=file)

		edit = Menu(pasek)
		edit.add_command(label="Undo")
		pasek.add_cascade(label="edit!", menu=edit)

		edit.add_command(label="tekst ;)", command=self.lol)

		canvv = Canvas(self.master, width=200, height=100)
		canvv.pack()

        canvv.create_rectangle(50, 25, 150, 75)

	def exit_program(self):
		exit()

	def lol(self):
		text = Label(self, text="Hey, good lookin!")
		text.pack()

root = Tk()
root.geometry(str(RozmiarX) + "x" + str(RozmiarY))

app = Window(root)
root.mainloop()