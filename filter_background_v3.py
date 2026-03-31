###################################################################################
#                                                                                 #
#  This software is developed in the frame of the IXPE collaboration              #
#  following the prescriptions reported in A. Di Marco et al., AJ 165 143 (2023)  #
#  for the IXPE data analysis                                                     #
#  Repository page: https://github.com/aledimarco/IXPE-background                 #
#  February, 28, 2025 released version 2.3 fixing issues due to astropy           #
#                                                                                 #
#  History: v0 May, 20, 2022 - preliminar version with energy independent cuts    #
#           v1 September, 16, 2022 - preliminar energy dependent cut              #
#           v2 March, 3, 2023 - final energy dependent cut published              #
#           v2.1 July, 20, 2023 - first public release                            #
#           v2.2 July, 30, 2024 - revision tag column                             #
#           v2.3 February, 28, 2025 - fixing issues due to astropy                #
#           v3.0 March, 31, 2026 - updated after DU2 anomaly with different       #
#                  approaches depending on DU and time                            #
#                                                                                 #
###################################################################################

import os
import sys
import argparse
from astropy import log
from astropy.logger import logging
from astropy.io import fits
from astropy import wcs
import numpy as np
from astropy.time import Time

KEYWORDS = ['TCTYP7', 'TCUNI7', 'TCRPX7', 'TCRVL7', 'TCDLT7',
            'TCTYP8', 'TCUNI8', 'TCRPX8', 'TCRVL8', 'TCDLT8',
            'TUNIT2', 'TUNIT5']

# rejection rules
def cut_pix(pi):
    ene=pi*0.04
    y=70+(ene)*30
    return y

def cut_fra(pi):
    ene=pi*0.04
    y=0.8*(1-np.exp(-(ene+0.25)/1.1))+ene*0.004
    return y

def cut_pix_du2(pi):
    ene=pi*0.04
    y=75+(ene)*41
    return y

def cut_fra_du2(pi):
    ene=pi*0.04
    y=0.71*(1-np.exp(-(ene+0.21)/1.5))+ene*0.0012
    return y

def rejection(path_lv2, path_lv1, output):
    #opening lv2
    data_lv2={}
    print(path_lv2,' opening')
    hdulist_input = fits.open(path_lv2)
    extension_names = list(map(lambda _ext: _ext.name, hdulist_input[1:]))
    events = hdulist_input['EVENTS'].data.T
    hdulist_input.info()
    _du = hdulist_input[1].header['DETNAM'] #storage of the DU name
    _obs_time = hdulist_input[1].header['DATE-OBS']
    _obs_time = Time(_obs_time, format="isot", scale="utc")
    print(_obs_time.mjd)
    status=False
    status=np.logical_and(_du=='DU2', _obs_time.mjd>60779)
    print(_du,' Date Observation:',_obs_time,' Status DU2 anomaly:',status)
    
    _run = hdulist_input[1].data
    for _key in ['TRG_ID','X', 'Y', 'Q', 'U', 'PI', 'TIME']:
        data_lv2[_key] = _run[_key]

    # stoage of unit information, mainly for wcs
    keywords_old = {}
    for _key in KEYWORDS:
        keywords_old[_key]=hdulist_input[1].header[_key]
    print(hdulist_input[1].columns)

    # opening lv1
    keys=['NUM_PIX','EVT_FRA','TRK_BORD','TIME']
    file_list = {}
    data_lv1 = {_key: np.array([]) for _key in keys}
    
    for _file in path_lv1:
        print(_file,' opening')
        hdul = fits.open(_file)
        extension = list(map(lambda _ext: _ext.name, hdulist_input[1:]))
        events_l1 = hdul['EVENTS'].data.T
        hdul.info()
        _run = hdul[1].data
        for _key in keys:
            data_lv1[_key]=np.append(data_lv1[_key], _run[_key])

    print('START {} and {}'.format(data_lv1['TIME'][0],data_lv2['TIME'][0]))
    print('STOP {} and {}'.format(data_lv1['TIME'][-1],data_lv2['TIME'][-1]))
    
    # build the dictionary to filter
    data_filt={}
    a,_mask,c=np.intersect1d(data_lv1['TIME'],data_lv2['TIME'],return_indices=True)
    for _key in data_lv1:
        data_filt[_key]=data_lv1[_key][_mask]
    a,_mask,c=np.intersect1d(data_lv2['TIME'],data_lv1['TIME'],return_indices=True)
    data_filt['PI']=data_lv2['PI'][_mask]
    data_filt['X']=data_lv2['X'][_mask]
    data_filt['Y']=data_lv2['Y'][_mask]

    # filtering
    if status==True:
        efra =  np.logical_and(data_filt['EVT_FRA']>cut_fra_du2(data_filt['PI']),data_filt['EVT_FRA']<0.9)
        numpix = np.logical_and(efra,data_filt['NUM_PIX']<cut_pix_du2(data_filt['PI']))
        trk_bord = np.logical_and(numpix,data_filt['TRK_BORD']<2)
    else:
        efra =  np.logical_and(data_filt['EVT_FRA']>cut_fra(data_filt['PI']),data_filt['EVT_FRA']<1.0)
        numpix = np.logical_and(efra,data_filt['NUM_PIX']<cut_pix(data_filt['PI']))
        trk_bord = np.logical_and(numpix,data_filt['TRK_BORD']<2)
        
    if output=='bkg':
        mask = np.where(np.logical_not(trk_bord))[0]
    else:
        mask = np.where(trk_bord)[0]

    times_filtered = data_filt['TIME'][mask]
    times_lv2 = events['TIME']
    intersect, lv2_idx, lv1_idx = np.intersect1d(times_lv2, times_filtered, return_indices=True)
    
    if output!='tag':
        dict_item = lambda _col: (_col, events[_col][lv2_idx])
        columns = list(map(lambda _col: _col.name, events.columns))
        dict_events = dict(map(dict_item, columns))

        #create new fits file
        _header_primary = hdulist_input[0].header
        primary_hdu = fits.PrimaryHDU(header=_header_primary)
        hdul_new = fits.HDUList([primary_hdu])
        for extname in extension_names:
            if extname == 'EVENTS':
                input_columns = hdulist_input[extname].columns
                columns = []
                for i, _input_column in enumerate(input_columns):
                    _name, _format = _input_column.name, _input_column.format
                    _array = dict_events[_name]
                    if _name == 'X':
                        columns.append(fits.Column(name=_name, array=_array, format=_format, unit='pixel', coord_type=keywords_old['TCTYP7'],
                                                   coord_unit=keywords_old['TCUNI7'], coord_ref_point=keywords_old['TCRPX7'],
                                                   coord_ref_value=keywords_old['TCRVL7'], coord_inc=keywords_old['TCDLT7']))
                    elif _name == 'Y':
                        columns.append(fits.Column(name=_name, array=_array, format=_format, unit='pixel', coord_type=keywords_old['TCTYP8'],
                                                   coord_unit=keywords_old['TCUNI8'], coord_ref_point=keywords_old['TCRPX8'],
                                                   coord_ref_value=keywords_old['TCRVL8'], coord_inc=keywords_old['TCDLT8']))
                    elif _name == 'TIME':
                        columns.append(fits.Column(name=_name, array=_array, format=_format, unit=keywords_old['TUNIT2']))
                    elif _name == 'PI':
                        columns.append(fits.Column(name=_name, array=_array, format=_format, unit=keywords_old['TUNIT5']))
                    else:
                        columns.append(fits.Column(name=_name, array=_array, format=_format))

                    
                _header_table = hdulist_input[extname].header
                table_hdu = fits.BinTableHDU.from_columns(columns, header=_header_table)

                hdul_new.append(table_hdu)
        
            else:
                hdul_new.append(hdulist_input[extname])

        hdul_new.writeto(path_lv2[:-5]+'_'+output+'.fits', overwrite=True)
        
        print(hdul_new[1].columns)

    else:
        #add masked column
        tag=np.zeros(len(data_lv2['TIME']))
        tag[lv2_idx]=1

        dict_item = lambda _col: (_col, events[_col])
        columns = list(map(lambda _col: _col.name, events.columns))
        dict_events = dict(map(dict_item, columns))

        #create new fits file
        _header_primary = hdulist_input[0].header
        primary_hdu = fits.PrimaryHDU(header=_header_primary)
        hdul_new = fits.HDUList([primary_hdu])
        for extname in extension_names:
            if extname == 'EVENTS':
                input_columns = hdulist_input[extname].columns
                columns = []
                for i, _input_column in enumerate(input_columns):
                    _name, _format = _input_column.name, _input_column.format
                    _array = dict_events[_name]
                    if _name == 'X':
                        columns.append(fits.Column(name=_name, array=_array, format=_format, unit='pixel', coord_type=keywords_old['TCTYP7'],
                                                   coord_unit=keywords_old['TCUNI7'], coord_ref_point=keywords_old['TCRPX7'],
                                                   coord_ref_value=keywords_old['TCRVL7'], coord_inc=keywords_old['TCDLT7']))
                    elif _name == 'Y':
                        columns.append(fits.Column(name=_name, array=_array, format=_format, unit='pixel', coord_type=keywords_old['TCTYP8'],
                                                   coord_unit=keywords_old['TCUNI8'], coord_ref_point=keywords_old['TCRPX8'],
                                                   coord_ref_value=keywords_old['TCRVL8'], coord_inc=keywords_old['TCDLT8']))
                    elif _name == 'TIME':
                        columns.append(fits.Column(name=_name, array=_array, format=_format, unit=keywords_old['TUNIT2']))
                    elif _name == 'PI':
                        columns.append(fits.Column(name=_name, array=_array, format=_format, unit=keywords_old['TUNIT5']))
                    else:
                        columns.append(fits.Column(name=_name, array=_array, format=_format))

                columns.append(fits.Column(name='BKG_TAG', array=tag, format='B'))
                _header_table = hdulist_input[extname].header
                table_hdu = fits.BinTableHDU.from_columns(columns, header=_header_table)

                hdul_new.append(table_hdu)
        
            else:
                hdul_new.append(hdulist_input[extname])

        hdul_new.writeto(path_lv2[:-5]+'_'+output+'.fits', overwrite=True)
        print(hdul_new[1].columns)
        
def main(args=None):
    
    parser = \
        argparse.ArgumentParser(description="Filter a lv2 file based on a lv1 cut")

    parser.add_argument("path_lv2", help="Input file lv2", type=str)
    parser.add_argument("path_lv1", help="Input file lv1", type=str, nargs='+')
    parser.add_argument("--output", default='rej',
                        help="Output file: rej (includes only source events), bkg (includes only rejected events) and tag (includes a tag column where events having 1 are src and 0 are bkg)", type=str)

    args = parser.parse_args(args)
    rejection(args.path_lv2, args.path_lv1, args.output)

if __name__ == '__main__':
    main(sys.argv[1:])
