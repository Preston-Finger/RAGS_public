"""
 Preston Mitchell Finger
 Minnesota State University, Mankato
 preston.finger@mnsu.edu, prestonfinger@gmail.com
"""

##########################################################################
##############################   Initialization  ##################
##########################################################################
# To activate from command line
# cl > ipython
# ipython > import RAGS as rags

# -----------------------------------Core Functionality ----------------------------------
import math
import os
#import sys
import pdb

import tkinter.ttk

import numpy as np
#from glob import glob
import time

#import importlib  ## utility in importlib.reload(insp) when debugging
#from distutils.sysconfig import *  # question-- what is this for?
import aplpy

#---------------------------------- Recent additions ---------------------------------------
#import timeit       # PMF used to see how long things take
import logging

# ---------------------------------------Tkinter for GUI----------------------------------
from tkinter import *
#from tkinter.ttk import *
from PIL import Image, ImageTk

from screeninfo import get_monitors

# ----------------------------------------Matplotlib--------------------------------------
import matplotlib
import matplotlib.pyplot as plt
#from matplotlib.figure import Figure
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)

import matplotlib.gridspec as gridspec

from scipy.optimize import curve_fit

#-----------------------------------------Astropy------------------------------------------
import astropy.io.fits as fits
from astropy.wcs import WCS as WCSa
#from astropy.coordinates import SkyCoord
#from astropy import units as u
from astropy.io import ascii
#from matplotlib.gridspec import GridSpec

################################## Unused Imports  #########################################################
#import _thread       # PMF trying this out.
#import threading

#import matplotlib.image as mpimg

#import shutil, tarfile, distutils, fileinput, scipy
#import matplotlib.transforms as mtransforms
#import matplotlib.patheffects as PathEffects
#from astropy.coordinates import match_coordinates_sky
#from astropy import wcs
#from astropy.nddata import utils
#from astropy.table import Table

#import importlib  ## utility in importlib.reload(insp)

###################################################################################################################
#################################  Backend Astrophysics  ######################################################
###################################################################################################################

## FIN in the line below needs to be fed after you read the catalog of objects you're going to survey

def plot_object(specid_, ra_, dec_, z_match1,Imag,Imag_error, passed):
    '''this is where the object gets plotted
    '''
    ra_ = float(ra_)
    dec_ = float(dec_)

    z_match1 = float(z_match1)
    h, w, px = Get_screen_size()

    fig = plt.figure(figsize=(w*px*1.25, h*px*0.74074))   # math to get it to fit in how we want
    #fig = plt.figure(figsize=(12, 6))  # Hard coded for MJR labtop
    #fig = plt.figure(figsize=(w * px * .85, h * px * 0.5))  # math to get it to fit in how we want
    matplotlib.rcParams['keymap.save'] ='x'

    # read the img
    sci = 'master_files/hlsp_candels_hst_wfc3_gs-tot_f160w_v1.0_drz.fits'
    #sci = 'master_files/hlsp_candels_hst_wfc3_gn-tot-60mas_f160w_v1.0_drz_nosip.fits'  # north nosip
    img = fits.open(sci)
    header = img[0].header
    data = img[0].data
    # read the positions

    ra, dec = ra_, dec_

    wcs_sol = WCSa(sci)
    y1, x1 = wcs_sol.wcs_world2pix(ra, dec, 0)
    w = WCSa(header)

    mnimg = np.mean(data[int(x1 - 10):int(x1 + 10), int(y1 - 10):int(y1 + 10)])
    # spec=fits.open('TESTSPEC/h_udf_wfc_id'+str(specid_)+'.fits')
    spec = fits.open('h_pears_total_south/h_pears_s_id' + str(specid_) + '.fits')

    wave, flux, ferror = spec[1].data['LAMBDA'], spec[1].data['FLUX'], spec[1].data['FERROR']

    headp = spec[1].header
    posange = float(headp['POSANG'].split('A')[1])
    spec.close()

    if not np.isnan(mnimg):
        rsline = (2175.0 * (1 + z_match1))

        if (np.min(wave) < rsline) & (np.max(wave) > rsline):

            newf = fits.PrimaryHDU()  # initialize a new object that will become a FITS file
            xlow, xhigh = int(x1 - 100), int(x1 + 100)
            ylow, yhigh = int(y1 - 100), int(y1 + 100)
            newf.data = img[0].data[xlow:xhigh, ylow:yhigh]  # this is the image to display,
            newf.header = img[0].header
            newf.header.update(w[xlow:xhigh, ylow:yhigh].to_header())
            newf.writeto('ID_' + specid_ + '_cutout_F160W.fits', overwrite=True)

            plt.gcf()

            gridspec_layout = gridspec.GridSpec(1, 3)
            gridspec_layout.update(hspace=0.0, wspace=0.0)

            gals = aplpy.FITSFigure('ID_' + specid_ + '_cutout_F160W.fits', figure=fig,
                                    subplot=[0.05, 0.35, 0.225, 0.275])
            galx = fits.open('ID_' + specid_ + '_cutout_F160W.fits')
            datag, headg = galx[0].data, galx[0].header
            galx.close()
            gmax, gmin = np.exp(-3) * np.max(datag), 0

            gals.axis_labels.set_xpad(100)
            gals.axis_labels.set_ypad(100)

            gals.show_grayscale(vmin=gmin, vmax=gmax)

            ax2 = fig.add_subplot(gridspec_layout[0, 1:])
            ax2.set(ylabel='Flux', xlabel='Wavelength(Angstroms)')  #todo stick y_label
            ax2.set_xlim(rsline - 1500, rsline + 1500)
            ax2.set_ylim(np.min(flux), np.max(flux))
            ax2.plot(wave, flux)  # creates the flux panel
            ax2.errorbar(wave, flux, yerr=ferror)
            ax2.text(0.5, 1, 'z = ' + str(z_match1), transform=ax2.transAxes)

            plt.axvline(rsline, color='r', linestyle='--')

            band_width = 150 * (1 + z_match1)
            plt.fill_between([rsline - band_width, rsline + band_width], np.min(flux), np.max(flux), alpha=0.5)

            grp = np.concatenate((np.where((wave > 1875 * (1+z_match1)) & (wave < 1975 * (1+z_match1))),
                                  np.where((wave > 2400 * (1+z_match1)) & (wave < 2500 * (1+z_match1))) ), axis=None)

            waver_, fluxr_ = wave[grp], flux[grp]

            popt, pcov = curve_fit(func, waver_, fluxr_)
            perr = np.sqrt(np.diag(pcov))
            print(waver_)
            print(pcov)
            print(perr)
            print(fluxr_)
            print_window(str(perr))
            ax2.plot(wave, func(wave, *popt), color='g')
            ax2.scatter(waver_, func(waver_, *popt), color='r')

            # PMF Hard coding for Heii 1640 Angstroms
            he2line = 1640. * (1+z_match1)
            smreg = np.where((wave > he2line-20) & (wave < he2line+20))

            if flux[smreg].size != 0:
                ax2.set_xlim(rsline - 2500, rsline + 1500)
                plt.axvline(he2line, color='c',
                        linestyle='--')  # or plt.axvline(rsline,'#008000') the # is a color hex value
                plt.fill_between([he2line - 20, he2line + 20], np.min(flux), np.max(flux), alpha=0.5)

            os.remove('ID_' + specid_ + '_cutout_F160W.fits')

    global canvas
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()

    ###############    TOOLBAR    ################ creating the Matplotlib toolbar
    frame = Frame(window)
    frame.grid(row=0,column=3)
    toolbar = NavigationToolbar2Tk(canvas, frame)
    toolbar.update()

    # placing the toolbar on the Tkinter window
    canvas.get_tk_widget().grid(row=2, column=0, sticky=W, columnspan=40)
    window.update()
    time.sleep(1)
    while (waiting):
        window.update()

    img.close()
    fig.clf()
    plt.clf()
    plt.close(fig)

def func(x, m,b):
    return m*x + b

def get_galaxies(file):
    ObjectID, RA, DEC, z_match, Imag, Imag_error = [], [], [], [], [], []
    with open(file, "r") as file:
        Lines = file.readlines()
        start = False
        c =  0
        for line in Lines:
            line = str.split((line))

            if starting_point[0] == "O" and line[0][0] == 'O':  # This is the starting point of new file with nothing in it
                start = True
                print("starting from beginning")

            elif starting_point == line[0]:
                print("We last left off at --> ", line[0])
                start = True
            elif start:
                ObjectID.append(line[0])
                RA.append(line[1])
                DEC.append(line[2])
                z_match.append(line[3])
                Imag.append(line[4])
                Imag_error.append(line[5])
                c += 1

        file.close()

    return ObjectID, RA, DEC, z_match, Imag, Imag_error, c

def inspect_objs():
    hide_widget(inspect_object)
    window.unbind('<Return>')

    window.update()

    label_directions = Label(text="Click 'a' for yes,  'n' for no  and 's' for maybe", font=("Times", 20))
    label_directions.grid(row=1, sticky=W)

    objid_unique, racks, dacks, zacks, Imag, Imag_error, c = get_galaxies("RAGS_input_dir/south_spectrum_matches.dat")

    speclist = "RAGS_input_dir/south_spectrum_matches.input"
    fins = ascii.read(speclist, names=['fnm'])

    donefile = 'done_%s' % user_name

    speclist_out = os.path.join('catalog_%s.dat' % user_name)
    if os.path.isfile(speclist_out):
        print('\nOutput file: \n  %s \nalready exists\n' % speclist_out)
        ask = input('Append? [Y/n] ')
        if ask.lower() == 'n':
            os.unlink(speclist_out)  # unlink method deletes the file
            os.unlink(donefile)
            objid_done = np.array([])
            # os.unlink(commentsfile)
            # starting over, no objects have been done
            os.unlink(donefile)
            objid_done = np.array([])
        else:
            objid_done = np.atleast_1d(np.genfromtxt(donefile, dtype=int))
    else:
        if os.path.exists(donefile):
            os.unlink(donefile)
        objid_done = np.array([], dtype=int)

    #### STEP 5:  Loop through objects ############
    #########################################################################

    for k in range(0, c):

        global waiting,label_object
        waiting=True

        progress = (len(fins['fnm']) - float(c) + k) / len(fins['fnm']) * 100.


        pb1['value'] = progress

        #next_obj = remaining_objects[0]
        current_object = objid_unique[k]

        tpl = [objid_unique[k], racks[k], dacks[k], zacks[k], Imag[k], Imag_error[k]]


        window.bind('a', lambda event: [write_to_file(tpl,'1')])   #yes
        window.bind('n', lambda event: [write_to_file(tpl,'0')])     #no
        window.bind('s', lambda event: [write_to_file(tpl,'2')])     #maybe

        window.bind('<Escape>', lambda event: warning_quit())

        hide_widget(label_object)
        label_object = Label(text='ObjectID: %i' % int(current_object), font=("Times",20, "bold"))
        label_object.grid(row=0, column=0, sticky=NE, columnspan=3)

        passval = 99

        classify(objid_unique[k],racks[k], dacks[k],  zacks[k],Imag[k],Imag_error[k] ,passval)
        window.update_idletasks()
    pb1['value'] = 100

    #clear_window()  # PMF  Clear the window of all widgets once done


    label = Label(master=window, font=('Bold', 50), text="You are done!\r Go eat some pizza")
    label.grid(row=3, column=3)

    msg = "Done with the main loop"
    print(msg)
    print_window(msg)

    window.unbind('a')    # PMF  Ok these buttons can go away now
    window.unbind('n')
    window.unbind('s')
    window.bind('<Return>', lambda event: window.quit())  # PMF close the window

def classify(ojacks, racks, dacks,  zacks,Imag,Imag_error, passval):
    print_window("...................Classifying.................")
    hide_widget(inspect_object)
    plot_object(ojacks,racks, dacks, zacks,Imag,Imag_error, passval)

    return

def get_remaining_objects(full_obj_list, objid_done):
    # array of flags
    tarr = np.array(full_obj_list)
    wdone = np.in1d(np.array(full_obj_list), objid_done)
    remaining = np.copy(tarr)
    mask = np.ones(tarr.shape, dtype=bool)
    mask[wdone] = False
    remaining = remaining[mask]
    return remaining

###################################################################################################################
#########################################  Utility Subroutines - Working On ########################################
###################################################################################################################
def delete_last_entry(): #todo figure this out
    window.unbind('a')
    window.unbind('n')
    window.unbind('b')

def help():
    #### NONE OF THIS IS ACTUALLY FUNCTIONAL
    msg = "Available Options:\n"
    msg += "\ta = accept source as absorber\n"
    msg += "\tac = accept with new shift applied\n"
    msg += "\ts = shift by the centroid by some amount\n"
    msg += "\tf = accept, but flag as possibly contaminated\n"
    msg += "\tr = reject as noise or irretrievably contaminated\n"
    msg += "\tq = quit and save"
    msg += "\th = I added this because why not"        #added line    also chang + - to just += for the msg's
    return msg

def print_vs_logging(): # PMF Maybe try out the logging import
    logging.debug("Debug info here")
    logging.info("just some info")
    logging.error("uh oh :(")

###################################################################################################################
#########################################  Utility Subroutines - Functional ########################################
###################################################################################################################
def Get_screen_size():
    h, w = math.floor(get_monitors()[0].height), math.floor(get_monitors()[0].width)
    px = 1 / plt.rcParams['figure.dpi']  # pixel in inches
    return int(h), int(w), px  #fixme remove the /1.5 as it is for troubleshooting

def print_window(msg):
    global text_box
    str(msg)
    text_box.insert(1.0, msg + "\n\n")   # PMF Changed this line
    return

def warning_quit():
    warning_window = Tk()
    warning_window.title('Warning')
    warning_window.geometry("200x200+500+500")
    window.bind('<Escape>', lambda event: quit_window)

    btn1 = Button(warning_window,text='Quit', command=quit_window)
    btn1.pack()
    btn2 = Button(warning_window,text='Continue', command=warning_window.destroy)
    btn2.pack()

def quit_window():
    global waiting
    waiting = False
    clear_window()
    window.quit()

def clear_window():                 # PMF This fucntion clears all widgets from window
    list = window.grid_slaves()
    for l in list:
        l.destroy()

def hide_widget(widget):  # This function is to remove widgets
    widget.grid_forget()
###################################################################################################################
#########################################  Utility Subroutines - File System ######################################
###################################################################################################################
def find_starting_point(file):
    global starting_point
    with open(file, "r") as file:
        last_line = file.readlines()[-1]
        last_line = str.split(last_line)
        file.close()
    starting_point = last_line[0]

def get_user():
    print_window("Ready to start inspecting objects?")

    global write_name
    window.unbind('<Return>')
    window.bind('<Return>', lambda event: inspect_objs())
    write_name = user_name.get()

    hide_widget(user_name_button)
    hide_widget((user_name))

    file_exists('users/south_Classification_%s.txt' %write_name)

    inspect_object.grid(row=0,column =0, pady=20)

def file_exists(file):
    file_exists = os.path.exists(file)

    if (file_exists is False):
        open(file, "a+")
        tpl = ["ObjectID", "RA", "DEC", "z_match", "Imag", "Imag_error"]
        write_to_file(tpl, 'Classification')
    else:
        print("You have file....lets keep going \n")

    find_starting_point(file)

def add_padding(tpl):
    line = ""
    for j in tpl:
        line += str(j).ljust(20, ' ')
    return line

def write_to_file(tpl, classification):
    window.unbind('a')  # PMF so we don't hit the key before the next galaxy spectrum
    window.unbind('n')
    window.unbind('s')

    tpl.append(classification)
    line = add_padding(tpl)

    msg = "writing -->  ", tpl[0], " to file as a ", classification,"\n"
    msg = ''.join(msg)
    print_window(msg)

    global waiting
    waiting = False

    canvas.get_tk_widget().grid_remove()
    window.update()

    user_file = open(('users/south_Classification_%s.txt' % write_name), "a+")
    user_file.write(str(line))
    user_file.write("\n")
    user_file.close()

###################################################################################################################
########################################## GUI Interface ##########################################################
###################################################################################################################

h, w, px = Get_screen_size()

window = Tk()      # PMF the main Tkinter window or just Tk window
window.title('Rapid Assessment of Galaxy Spectra')     # PMF setting the title of the Tk window

window.geometry(str(w) +'x'+ str(h) + '+' + "0+0")

load = Image.open("images/milky-way.jpg")
bg = ImageTk.PhotoImage(load)

background = Label(window,image=bg)
background.place(x=0,y=0, relwidth=1, relheight=1)

waiting = True   # PMF Global variable so we wait until galaxy is classified

#   PMF These are the Tk various widgets.  Buttons, labels, entrys, etc.
#   They are all set here so every function can use them.

window.bind('<Escape>', lambda event: warning_quit)

label_position = Label(master=window)
label_object = Label(master=window)
label_msg = Label(master=window)

user_name_button = Button(window, text="Enter your initials all caps", command=get_user, font=("Bold", 20))
user_name_button.grid(row=0, column=0, pady=1)
window.bind("<Return>", lambda event: get_user())

user_name = Entry(master=window)
user_name.grid(row=1, column=0, pady=1)
user_name.focus_set()                       # PMF this starts the cursor in the Entry box at start.

inspect_object = Button(master=window, command=inspect_objs, text="Inspect Object", font=("Bold", 50))

text_box = Text(window, width = int(w/19.2), height = int(h/27))
text_box.grid(row=3, column=0, sticky=W)

scroll_y = Scrollbar(window, orient = "vertical", command=text_box.yview)
scroll_y.grid(row=3, column=0, sticky=E)

pb1 = tkinter.ttk.Progressbar(window, orient=HORIZONTAL ,  length=400, mode='determinate', style='TProgressbar')
pb1.grid(row=3, column=1)

fig = plt.figure(figsize=(w*px*1.25, h*px*0.74074))   # math to get it to fit in how we want
#fig = plt.figure(figsize=(12, 6))  # Hard coded for MJR labtop
global canvas
canvas = FigureCanvasTkAgg(fig, master=window)
canvas.draw()
canvas.get_tk_widget().grid(row =2, column=0,sticky=W, columnspan=40)


print_window("Enter in your 3 initials")

# run the gui
window.mainloop()