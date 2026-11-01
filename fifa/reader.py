#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 2026
@author: lucindawang
"""

# Standard library imports
import re
import os
import astropy 
import openpyxl

# Basic package imports
import numpy as np 
from matplotlib import lines
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
np.trapz = trapezoid 

"""
Welcome to FIFA!!
Here we provide two main scripts to read force-photometry data from ZTF and ATLAS respectively.
After you earn your own data file from each of the survey. Please insert the path of your photometry file in the following line:
"""

class SN_reader_ztf():
    """
    Parameters
    ----------
    dpath : str
        The path to the ZTF data file (e.g., self.dpath = f'/the/way/to/your/ztffps_file.txt')
    n : int
        The line number to read from the ZTF data file (0-indexed).
    """
    def __init__(self, dpath):
        self.dpath = dpath
    """
    1. The script used to read the line on by one. Firstly, we replace the space with comma.
    2. Find lines containing target_word (skip other # lines). 
    3. Skip all lines containing #, but keep target_line.
    """
    def read_line(self, n:int): 
        with open(self.dpath, 'r') as f:
            file = f.readlines()
            new_data = file[n].replace(" ", ",").split(",") 
            clean_data = [str(x) for x in new_data if x != '']
        return clean_data
    
    def read_index(self): 
        with open(self.dpath, "r", encoding = "utf-8") as file:
            lines = file.readlines()
            target_word = '#'
            positions   = np.array([])
        for line_number, line in enumerate(lines, start = 0): 
            if target_word in line: 
                pass
            else:
                positions = np.append(positions, line_number)
        return positions

    def read_data_index(self): 
        with open(self.dpath, "r", encoding = "utf-8") as file:
            lines = file.readlines()
            target_word = '#'
            positions = np.array([])

        for line_number, line in enumerate(lines, start = 0): #originally start = 0
            start = 0
            while True:
                index = line.find(target_word, start)
                if index == -1:
                    break
                positions = np.append(positions, (line_number, index+1))
                start     = index + len(target_word)
        last_1 = int(positions[-1][0])
        last_2 = int(positions[-2][0])
        return last_2, last_1

    def earn_data(self, key=None): 
        SNT = 3 #signal-to-noise ratio
        reader      = SN_reader_ztf(self.dpath)
        index       = reader.read_index() #列出所有＃的位置
        data_index  = reader.read_data_index() #只列出最後兩個＃的位置
        index_array = reader.read_line(index[0]) 
        end_num     = data_index[1]
        start_num   = data_index[0]

        col = {name: i for i, name in enumerate(index_array)}

        #Read the row data once, avoid repeat calling
        row_line = [reader.read_line(i) for i in range(start_num + 1, end_num)]
        sn_mjd, sn_filter, sn_zp, sn_flux, sn_fluxerr = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
        """
        Quality Filtering
        -----------------
        These were proven to be useful following analyses of large numbers of Type-Ia SN lightcurves (e.g., Yao et al. 2019).
        1. remove forcediffimflux_box == 'null'
        2. remove info_box >= 33554432 (=2^(25))
        3. remove sci_box >= 25 (DN)
        4. remove seeing_box >= 4 (arcsec)
        5. remove SNT_box < 3
        """
        for data in row_line:

            if data[col["forcediffimflux"]] == 'null' or data[col["forcediffimfluxunc"]] == 'null':
                continue

            f   = float(data[col["forcediffimflux"]])
            ferr= float(data[col["forcediffimfluxunc"]])
            
            if ferr <= 0 or (f/ferr) <= 3:
                continue
            
            info   = float(data[col["infobitssci"]])
            sci    = float(data[col["scisigpix"]])
            seeing = float(data[col["sciinpseeing"]])
            if info >= 33554432 and sci >= 25 and seeing >= 4:
                continue
            
            band = data[col["filter"]].split('_')[1]
            sn_filter = np.append(sn_filter, f'ztf_{band}')
            sn_mjd    = np.append(sn_mjd, float(data[col["jd"]]) - 2400000.5)
            sn_zp     = np.append(sn_zp, float(data[col["zpdiff"]]))
            sn_flux   = np.append(sn_flux, f)
            sn_fluxerr= np.append(sn_fluxerr, ferr)
        """ 
        For the details on why the SNT=3, the related paper is given by this:
        http://web.ipac.caltech.edu/staff/fmasci/home/mystats/UpperLimits_FM2011.pdf
        """
        mask  = ((sn_flux/sn_fluxerr) > SNT)
        mag   = sn_zp[mask] - 2.5*np.log10(sn_flux[mask])
        sigma = 1.0857*(sn_fluxerr[mask]/sn_flux[mask])

        data_dict = {
            'mjd': sn_mjd,
            'filter': sn_filter,
            'flux': sn_flux,
            'flux_unc': sn_fluxerr,
            'zp': sn_zp,
            'mag':mag,
            'sigma':sigma
        }
        if key is None:
            return data_dict
        try:
            return data_dict[key.lower()]
        except KeyError:
            raise KeyError(f"Invalid key '{key}'. Available keys: {list(data_dict)}")
    
    
    def run_figure(self, sn_name:str):
        reader = SN_reader_ztf(self.dpath)
        time   = reader.earn_data()['mjd']
        band   = reader.earn_data()['filter']
        mag    = reader.earn_data()['mag']
        err    = reader.earn_data()['sigma']

        fmt_list = ['o', 's', '^', 'D', 'v', 'p', '*', 'h', 'H', 'X']
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
        mpl.rcParams['font.family'] = "monospace" #"Times New Roman"
        mpl.rcParams['lines.linewidth'] = 10.0
        for i, b in enumerate(np.unique(band)):
            plt.errorbar(time[band == b], mag[band == b], label = f'{b[4]}-band', yerr = err[band == b], fmt = fmt_list[i % len(fmt_list)], ecolor = 'grey', color = colors[i % len(colors)], elinewidth = 1, capsize = 3)
            plt.ylim(np.max(mag[band == b])+1, np.min(mag[band == b])-1)
            plt.legend()
            plt.ylabel('Apparent Magnitude')
            plt.xlabel('MJD')
            folder_path = os.path.dirname(self.dpath)
            new_folder = os.path.join(folder_path, sn_name)
            os.makedirs(new_folder, exist_ok=True)
            plt.title(f'griz Light curve of SN{sn_name}'.format(sn_name))
            plt.savefig(f'{new_folder}/SN{sn_name}_ztffps.png', dpi=600, facecolor='w', edgecolor='w', orientation='portrait')
        return plt.show()

"""
Now, we have to read the data from ATLAS FPS file, please insert the path to your ATLAS FPS file in the following line:
"""
class SN_reader_atlas():
    """
    Parameters
    ----------
    dpath : str
        The path to the ATLAS data file (e.g., self.dpath = f'/the/way/to/your/atlasfps_file.txt')
    n : int
        The line number to read from the ATLAS data file (0-indexed).
    """
    def __init__(self, dpath):
        self.dpath = dpath

    def read_line(self, n:int):
        f = open(self.dpath, 'r')
        file = f.readlines()
        if n == 0:
            data = file[:][n][1:]
            data2 = data.replace(" ", ",").split(",")
            data2 = [item for item in data2 if item.strip() != '']
        else:
            data = file[:][n][:]
            data2 = data.replace(" ", ",").split(",")
            data2 = [item for item in data2 if item.strip() != '']
        f.close()
        return data2
    
    def earn_data(self, key=None):
  
        index_value  = SN_reader_atlas(self.dpath).read_line(0)
        MJD_index    = index_value.index('##MJD') #set position for data
        uJy_index    = index_value.index('uJy')
        duJy_index   = index_value.index('duJy')
        ftype_index  = index_value.index('F')
        err_index    = index_value.index('err')
        x_index      = index_value.index('x')
        y_index      = index_value.index('y')
        maj_index    = index_value.index('maj')
        min_index    = index_value.index('min')
        apfit_index  = index_value.index('apfit')
        mag5sig_index= index_value.index('mag5sig')
        Sky_index    = index_value.index('Sky')

        # 1. Read a file once, deal with the whole row data
        with open(self.dpath, 'r') as f:
            lines = f.readlines()

        # Deal with the first row and split the following strip
        header = lines[0].lstrip('#').split()
        
        # Split the whole data
        data_rows = [line.split() for line in lines[1:] if line.strip()]

        # Trensfer into NumPy array, get the column data once
        data_matrix = np.array(data_rows)

        ftype_box   = data_matrix[:, ftype_index]
        err_box     = data_matrix[:, err_index].astype(float).astype(int)
        x_box       = data_matrix[:, x_index].astype(float)
        y_box       = data_matrix[:, y_index].astype(float)
        mag5sig_box = data_matrix[:, mag5sig_index].astype(float)
        MJD_box     = data_matrix[:, MJD_index].astype(float) #Firstly, make the float data
        maj_box     = data_matrix[:, maj_index].astype(float)
        min_box     = data_matrix[:, min_index].astype(float)
        uJy_box     = data_matrix[:, uJy_index].astype(float)
        duJy_box    = data_matrix[:, duJy_index].astype(float).astype(int)
        apfit_box   = data_matrix[:, apfit_index].astype(float)
        Sky_box     = data_matrix[:, Sky_index].astype(float)

        #criteria
        mask  = ((duJy_box<10000)&(err_box==0)& 
                (x_box>100)&(x_box<10460)&
                (y_box>100)&(y_box<10460)&
                (maj_box<5)&(maj_box>1.6)&
                (min_box<5)&(min_box>1.6)&
                (apfit_box>-1)&(apfit_box<-0.1)&
                (mag5sig_box>17)&(Sky_box>17))

        date_array   = np.array(MJD_box[mask])
        duJy_array   = np.array(duJy_box[mask])
        flux_array   = np.array(uJy_box[mask])
        type_array   = np.array(ftype_box[mask])
        unique_day   = np.unique((date_array).astype(int)) #beginning with the average part
        avg_type_box = np.zeros(len(unique_day), dtype = str)
        avg_m_box    = np.zeros(len(unique_day), dtype = float) #make average flux value for ezch day
        avg_sigma_box= np.zeros(len(unique_day), dtype = float)
        avg_err_box  = np.zeros(len(unique_day), dtype = float) 
        avg_day_box  = np.zeros(len(unique_day), dtype = float) 
        for i in range(len(unique_day)): 
            day_loc = np.where((date_array).astype(int) == unique_day[i]) #calculate how many times 
            n = len(day_loc[0])
            final_date = date_array[day_loc]
            final_type = type_array[day_loc]
            final_flux = flux_array[day_loc] #how many flux per day
            sigma  = duJy_array[day_loc] #how many errors per day
            weight = 1/(sigma*sigma)
            weight_avg  = np.sum(weight*final_flux)/np.sum(weight) #Weighted average
            weight_std  = 1/np.sqrt(np.sum(weight)) #Weighted standard deviation
            weight_date = np.sum(final_date)/n #Average date

            avg_m_box[i]     = weight_avg
            avg_err_box[i]   = weight_std
            avg_day_box[i]   = weight_date
            avg_sigma_box[i] = np.abs((2.5 / np.log(10)) * (weight_std / weight_avg)) #After I fight with my Gemini, I got this.
            avg_type_box[i]  = final_type[0]

        sigma_3 = 23.9 - 2.5*np.log10(3*avg_sigma_box)
        o_loc = np.where((avg_type_box == 'o') & (avg_m_box > 3*sigma_3)) #make criteria for o-band data
        c_loc = np.where((avg_type_box == 'c') & (avg_m_box > 3*sigma_3))

        odate = avg_day_box[o_loc]
        cdate = avg_day_box[c_loc]
        omag  = 23.9 - 2.5*np.log10(avg_m_box[o_loc])
        cmag  = 23.9 - 2.5*np.log10(avg_m_box[c_loc])
        oflux = avg_m_box[o_loc]
        cfluc = avg_m_box[c_loc]
        oflux_err = avg_err_box[o_loc]
        cflux_err = avg_err_box[c_loc]
        omag_err  = avg_sigma_box[o_loc]
        cmag_err  = avg_sigma_box[c_loc]

        data_dict = {
            'odate': odate,
            'cdate':cdate,
            'omag':omag,
            'cmag':cmag,
            'omag_err':omag_err,
            'cmag_err':cmag_err
        }
        if key:
            return data_dict[key.lower()]

        return data_dict
    
    def run_figure(self, sn_name):
        reader= SN_reader_atlas(self.dpath)
        odate = reader.earn_data()['odate']
        cdate = reader.earn_data()['cdate']
        omag  = reader.earn_data()['omag']
        cmag  = reader.earn_data()['cmag']
        omag_err = reader.earn_data()['omag_err']
        cmag_err = reader.earn_data()['cmag_err']
        mpl.rcParams['font.family'] = "monospace" #"Times New Roman"
        mpl.rcParams['lines.linewidth'] = 10.0
        with plt.style.context(['science', 'high-vis']):
            plt.errorbar(odate, omag, label = 'o-band', yerr = omag_err, fmt = 'o', ecolor = 'grey', color = 'b', elinewidth = 2, capsize = 4)
            #The part is for c-band
            plt.errorbar(cdate, cmag, label = 'c-band', yerr = cmag_err, fmt = 'o', ecolor = 'grey', color = 'g', elinewidth = 2, capsize = 4)
            plt.ylim(np.max(omag)+2,np.min(omag)-2)
            plt.grid()
            plt.xlabel('MJD')
            plt.ylabel('Apparent Magnitude')
            plt.title('oc Light curve of '+ sn_name)
            
            folder_path = os.path.dirname(self.dpath) #Find the position that the file located
            new_folder = os.path.join(folder_path, sn_name) #Create the new folder at the position
            os.mkdir(new_folder)

            plt.savefig(f'{new_folder}/SN{sn_name}_atlasfps.png', dpi=600, facecolor='w', edgecolor='w', orientation='portrait')
            plt.legend(loc = 'upper right')
        return plt.show()