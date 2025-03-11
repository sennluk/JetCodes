# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 17:12:41 2021

@author: egiov
"""

#import tkinter
import eqdsk
from tkinter import ttk, Tk, Button, N, W, E, S, filedialog

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

import numpy as np

class EqdskTk:
    def __init__(self, root):
        root.title("EqdskTk")
        mainframe = ttk.Frame(root, padding="3 3 12 12")
        mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
        self.canvas = self.create_canvas(mainframe)
        
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        
        mainframe.columnconfigure(0, weight=1)
        mainframe.rowconfigure(0, weight=1)
        
        button_frame = ttk.Frame(mainframe)
        button_frame.grid(row=2, sticky=(N,S,E,W))
        button_frame.columnconfigure(0, weight=1)
        button_frame.columnconfigure(1, weight=1)
        
        button_load = Button(master=button_frame, text='load',
                             command=self.load_eqdsk)
        button_load.grid(row=0, column=0, sticky=(W,E))
        
        button = Button(master=button_frame, text="Quit", command=root.quit)
        button.grid(row=0, column=1, sticky=(E,))
        
    def load_eqdsk(self):
        filename = filedialog.askopenfilename()
        if filename is None:
            return
        
        try:
            w = eqdsk.eqdsk_read(filename)
        except ValueError:
            print("Error")
            return
        print("lettura")
        self.fig.suptitle(w.case)
        first = True
        eqdsk.ploteqdsk_psi(w, self.ax[0], first=first)
        eqdsk.ploteqdsk_lin(w, self.ax[1], first=first)
        self.canvas.draw()
        
    def create_canvas(self, root):
        self.fig, self.ax = eqdsk.create_eqdsk_fig()
        #t = np.arange(0, 3, .01)
        #fig.add_subplot().plot(t, 2 * np.sin(2 * np.pi * t))

        canvas = FigureCanvasTkAgg(self.fig, master=root)  # A tk.DrawingArea.
        canvas.get_tk_widget().grid(row=0, sticky=(N,W,E,S))

        # pack_toolbar=False will make it easier to use a layout manager later on.
        toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False)
        toolbar.update()


        canvas.mpl_connect(
            "key_press_event", lambda event: print(f"you pressed {event.key}"))
        canvas.mpl_connect("key_press_event", key_press_handler)

        toolbar.grid(row=1, sticky=(W,E))
        return canvas

def main():
    root = Tk()
    EqdskTk(root)
    root.mainloop()
    
if __name__ == "__main__":
    main()

