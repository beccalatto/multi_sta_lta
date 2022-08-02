# -*- coding: utf-8 -*-
"""
__update__ = "Tues July 15, 2022"
__author__ = "Rebecca Latto"
__version__ = "Python 3.7"
__application__ = "Spyder"
__email__ = "rebecca.latto@utas.edu.au"
"""

# In[]: General imports

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from obspy import Stream, UTCDateTime
import matplotlib
import math 
import matplotlib.backends.backend_pdf
from matplotlib.ticker import MultipleLocator
import warnings
import os.path
from time import strftime
from time import gmtime

# In[]: seismic_attributes import

import seismic_attributes as sa

# In[Event confidence routine]: 
# Developed by Rebecca Latto and Ross J. Turner (<ross.turner@utas.edu.au>)
# Assign each station a list of its 5 nearest neighbors
# Check if any stations that detected an event are not 
# in nearest neighbor groups (i.e. are outlier, so likely spurious detections)
    #Then, remove that seismometer iteratively
        #If the remaining list of seismometers falls below the coincidence sum (here, 3), 
        # flag that event low confidence
    #Else, flag that event high confidence

from obspy.geodetics import locations2degrees

def confidence(event_stations,coincidence_sum,station_dict):
    
    station_nearest_neighbors = {}
    for key in station_dict:
        #print(key)
        WIS_lat = station_dict[key][0]
        WIS_lon = station_dict[key][1]  
        nearest_neighbor_list = []
        station_list = []
        for other_keys in range(len(station_dict)):
            nearest_neighbor_list.append(locations2degrees(list(station_dict.values())[other_keys][0],\
                                                           list(station_dict.values())[other_keys][1],\
                                                               WIS_lat,WIS_lon))
            station_list.append(list(station_dict.items())[other_keys][0])
            nearest_neighbors_list_ordered = [x for _, x in sorted(zip(nearest_neighbor_list, station_list), key=lambda pair: pair[0])]
        station_nearest_neighbors[key] = ((nearest_neighbors_list_ordered[1:6]))
    
    confidence_marker = _confidence(event_stations,4,station_nearest_neighbors)
    return confidence_marker

def _confidence(event_stations,coincidence_sum,station_nearest_neighbors):
    
    confidence_list = []
    for station_idx in range(len(event_stations)):
        temp_station = event_stations[station_idx]
        temp_station_NN = station_nearest_neighbors[temp_station]
        confidence_list.append(len([w for w in temp_station_NN if w not in event_stations]))
    
    # Does the confidence list contain any stations that are outliers
    # Outliers are stations that are not in any nearest neighbor groups
    # , so the confidence_list[value] would equal the length of the station list
    confidence_list_binary = []
    for i in range(len(confidence_list)):
        if confidence_list[i] >= len(event_stations):
            confidence_list_binary.append(0)
        else:
            confidence_list_binary.append(1)
    
    if sum(confidence_list_binary) >= coincidence_sum:
        confidence_marker = 'High'
    else: 
        confidence_marker = 'Low'
    print(confidence_marker)
    
    return confidence_marker

# In[]:
def myAnalystPlots(event_catalogue,stream_dict=None,station_dict=None,\
                   filters=True,highpass=100,lowpass=0.001,trace_order=None,\
                   component=None,vertical_scaling_range=1.5e3,\
                   arrival_buffer=60,total_duration=60*10,event_color='b',\
                   event_alpha=0.4,folder=None, max_number_PDF_slides = 20,\
                   confidence_marker=False,line_width=0.5,sampling_interval=4,\
                   tick_lines_every_x_secs=5,interactive=False,LabelList=None,\
                       label=None):
    
    #### event_catalogue: Panda Dataframe produced using seismic_attributes.py with  
    ####                  a reference catalogue and a trace catalogue
    #### stream_dict: Dictionary of stream attributes including network, station,
    ####              location, channel, start_time, end_time
    ####
    #### station_dict: Dictionary of stations, lat and lon locations
    ####
    #### filter: True or false if filter is to be applied
    ####
    #### highpass: User can define a highpass filter. Default is 100 Hz.
    ####
    #### lowpass: User can define a lowpass filter. Default is 0.001 Hz.
    ####
    #### trace_order: Traces can be ordered in the following three ways
    ####                1. None: indicates no ordering, and the traces will be 
    ####                    thus ordered in the order they appear in the input_stream
    ####                2. 'arrival_time': indicates order of traces depends on 
    ####                    the arrival time per event
    ####                3. [list]: List as input indicates that the order of traces
    ####                   is determined by the input list. Suggested, list is ordered
    ####                   consistently based on location of stations.
    ####
    #### component: Components can be plotted in the following ways:
    ####            - 'Z' : Just plot the zed
    ####            - 'E' : Just plot the east
    ####            - 'N' : Just plot the north
    ####            - 'ALL': Plot all three components, stacked
    ####
    #### vertical_scaling_range: (int) Default corresponds to ObsPy default (1.5e3). 
    ####                         Need to adjust depending on pre-processing filters or 
    ####                         detrending. A single value suggested for plot-to-plot
    ####                         consistency. However, argument None can be used to 
    ####                         invoke a variable y-limit choice.
    ####
    #### arrival_buffer: (int) Number of seconds before plotting the first arrival of
    ####                 an event. Default is one minute (60 seconds)
    ####
    #### total_duration: (int) Default is to have each plot total 10 minutes on the x-axis.
    ####                 , meaning 9 minutes following the default arrival buffer 
    ####                 regardless of event length (longer or shorter than 10 min)
    ####
    #### event_color: Events are highlighted on plots. Default is blue 'b' using
    ####                Python color codes
    ####
    #### event_alpha: Events are highlighted on plots. Default is 0.4 for a very
    ####                transparent highlight. Can be lighter (to 0.0) or darker (1.0).
    #### 
    #### folder: Folder in which user saves the .pdf outputs. Default is current
    ####         working directory.
    ####
    #### max_number_PDF_slides: (int) The pdf outputs are partitioned because of 
    ####                        performance capacity of Adobe PDF. Default is 20 
    ####                        based on author experience of Adobe PDF.
    ####                    
    #### confidence_marker: (boolean, True or False)
    ####            If False, nothing is changed. If True, calls the confidence()
    ####            routine which marks events by likelihood of real event 
    ####            based on a location metric. Plot will have an annotation saying
    ####            Low or High confidence depending on the routine output.
    ####
    #### line_width: (float) Default is 0.5. Python default is 1.5.
    ####
    #### sampling_interval: (int) Default is every 4th sample for plotting.    
    ####                    A higher number creates a finer line and a lower 
    ####                    number (minimum=1) creates a thicker line. Needs to
    ####                    be below the Nyquist frequency
    ####    
    #### tick_lines_every_x_secs: (int) Default is to include a tick marker every
    ####                            5 seconds.
    ####
    #### interactive: (boolean, True or False) Interactive mode has two options:
    ####            If False, the script will run as is and will produce .pdfs
    ####            and close iteratively produced plots.
    ####            If True, the script will pause after producing each plot to 
    ####            allow for analyst to zoom, move around, etc interactively.
    ####            When done, the user can save or just close the plot and the
    ####            next one will appear, and so on.
    ####
    #### LabelList: (list or None) If user has a list of event labels that corresponds
    ####            to reference event catalogue (event_catalogue[0]), then it can
    ####            can be included if user wants to sort events in MyAnalystPlots 
    ####            by label groups. Default is None.
    ####
    #### label: (string or None) If a LabelList is included as input, then the user
    ####         has the option to indicate which label a group should be sorted by.
    ####         For example, if label is 'Teleseism', the routine will only pull
    ####         the 'Teleseism' events from event_catalogue[0] using LabelList.
    ####         Default is None.
    
    warnings.filterwarnings("ignore")
    
    # Hard coded coincidence sum for Event detection for cryoseismology
    coincidence_sum = 3

    station_list = []
    for key in range(len(station_dict)):
        station_list.append(list(station_dict.items())[key][0])

    # Check for label filter (i.e. compare number of events with label to max_number_PDF_slides)
    try:
        count = list(LabelList).count(label)
    except: # Would reach except if LabelList or label is None
    #   then, user will produce plots for the whole catalogue.
        count = len(event_catalogue[0])
        
    extensions = math.ceil(count/max_number_PDF_slides)
    
    print(extensions)
    ext_counter = -1*max_number_PDF_slides
    
    print('Creating set(s) of pdfs based on the max_number_PDF_slides argument')
    
    event_times = event_catalogue[0]['ref_time']
    for iext in range(extensions):
        print('Creating pdf set %s' %(iext+1))
        #print(ext_counter)
        ext_counter += max_number_PDF_slides
        
        try: 
            ev_idxs = [i for i, x in enumerate(list(LabelList)) if x == label]
        except:
            ev_idxs = np.arange(len(event_catalogue[0]))

        # Is user creating just one set of .pdfs
        if extensions == 1:
            extension = ''
            EV_IDXS = ev_idxs
        # or, several sets of .pdfs
        else:
            extension = '_set%s' %iext
            EV_IDXS = ev_idxs[ext_counter:ext_counter+max_number_PDF_slides]
            print(EV_IDXS)
    
        # t1
        start_year = UTCDateTime(event_times[EV_IDXS[0]]).year
        start_mo = UTCDateTime(event_times[EV_IDXS[0]]).month
        start_day = UTCDateTime(event_times[EV_IDXS[0]]).day
        start_t = UTCDateTime(start_year,start_mo,start_day,0,0)
        # t2
        end_year = UTCDateTime(event_times[EV_IDXS[len(EV_IDXS)-1]]).year
        end_mo = UTCDateTime(event_times[EV_IDXS[len(EV_IDXS)-1]]).month
        end_day = UTCDateTime(event_times[EV_IDXS[len(EV_IDXS)-1]]).day
        end_t = UTCDateTime(end_year,end_mo,end_day,0,0)
        print(start_t,end_t)
        
        if interactive is False:
            version_extension = ''
            version_counter = 0
            file_exists = True
            while file_exists == True:
                try:
                    filename = 'MyAnalystPlot_%s%s%s.pdf' %(label,extension,version_extension)
                except: 
                    filename = 'MyAnalystPlot_%s%s.pdf' %(extension,version_extension)
                
                file_exists = os.path.exists(filename)   
                version_extension = '_v%s' %(version_counter)
                version_counter += 1
            
            # CREATES THE PDF
            pdf1 = matplotlib.backends.backend_pdf.PdfPages(filename)
        
        #plt.ion() indicates interactive on for matplotlib
        plt.ion()
        # Loop through each event
        print('Looping through each event in EV_IDXS')
        loop_idx = 0
        print(EV_IDXS)
        if iext == 0:
            continue
        for iev in EV_IDXS:
            print(event_times[iev])
            print(iev)
            new_year = UTCDateTime(event_times[iev]).year
            new_mo = UTCDateTime(event_times[iev]).month
            new_day = UTCDateTime(event_times[iev]).day
            new_t = UTCDateTime(new_year,new_mo,new_day,0,0)
            print(new_t)
                            
            print('Loading in stream objects while the date of iev is not the same as previous')
            if loop_idx == 0: # First one will be equal to start_t
                input_stream = Stream()
                for i in range(0,len(station_list)):
                    st = Stream()
                    station = station_list[i]
                    # Max loads one day at a time
                    st = sa.get_waveforms(stream_dict['network'], station, \
                                          stream_dict['location'], \
                                          stream_dict['channel'], \
                                          start_t, start_t+86400, \
                                          event_buffer=3600, \
                                          waveform_name='waveforms',\
                                          station_name='stations',\
                                          download=False)    
                    input_stream += st.copy().detrend(type='demean')
                
                if filters == True:
                    print('Filtering according to default or user input')
                    input_stream = input_stream.filter("bandpass", freqmin=lowpass,\
                                        freqmax=highpass,corners=4, zerophase=True)
                    
            elif new_t != start_t: # If the new event time does not equal start_t, change start_t
                print("new start_t is %s" %new_t)
                start_t = new_t
                input_stream = Stream()
                for i in range(0,len(station_list)):
                    st = Stream()
                    station = station_list[i]
                    # Max loads one day at a time
                    st = sa.get_waveforms(stream_dict['network'], station, \
                                          stream_dict['location'], \
                                          stream_dict['channel'], \
                                          start_t, start_t+86400, \
                                          event_buffer=3600, \
                                          waveform_name='waveforms',\
                                          station_name='stations',\
                                          download=False)    
                    input_stream += frequency_normalisation(st.copy().detrend(type='demean'))  
                
                if filter == True:
                    print('Filtering according to default or user input')
                    input_stream = input_stream.filter("bandpass", freqmin=lowpass,\
                                        freqmax=highpass,corners=4, zerophase=True)
            elif new_t == start_t:
                input_stream = input_stream # Already loaded
            
            loop_idx += 1
            
            ev_id = event_catalogue[0]['event_id'][iev]
            ref_dur = event_catalogue[0]['ref_duration'][iev]
            ev_id_idx = event_catalogue[1]['event_id'][event_catalogue[1]['event_id'] ==ev_id].index
            
            stations_that_det_event = list(event_catalogue[1]['station'][ev_id_idx]) 
             
            ev_t1 = event_catalogue[1]['time'][ev_id_idx]  
            ev_dur = event_catalogue[1]['duration'][ev_id_idx]
            earliest_arrival_idx = np.argmin(ev_t1)
    
            earliest_arrival = UTCDateTime(ev_t1[ev_id_idx[earliest_arrival_idx]])
            LHS_t1 = earliest_arrival - arrival_buffer
            RHS_t2 = earliest_arrival + total_duration - arrival_buffer
                
            # trace order
            if trace_order == None: # defaulted
                reordered_E = np.arange(0,len(input_stream),3)
                reordered_N = np.arange(1,len(input_stream),3)
                reordered_Z = np.arange(2,len(input_stream),3)
                reordered_station_list = station_list
                trace_order_sh = 'default'
                
            elif type(trace_order) == list: # trace order defined
                E = np.arange(0,len(input_stream),3)
                d = {v: k for k, v in enumerate(trace_order)}
                trace_order_mapping = list(map(d.__getitem__, station_list))
                reordered_E = [E[ii] for ii in trace_order_mapping]
                reordered_N = [ii+1 for ii in reordered_E]
                reordered_Z = [ii+2 for ii in reordered_E]
                reordered_station_list = [station_list[ii] for ii in trace_order_mapping]
                trace_order_sh = 'station location'
                
            elif type(trace_order) == str: # by arrival time
                E = np.arange(0,len(input_stream),3)
                stations_arrival = []
                TIMEs = event_catalogue[1]['time'][ev_id_idx]
                ev_id_idx_copy = list(ev_id_idx)
                for ii in range(len(event_catalogue[1]['station'][ev_id_idx])):
                    find_earliest_per_ii = np.argmin(TIMEs)
                    stations_arrival.append(event_catalogue[1]['station'][ev_id_idx_copy[find_earliest_per_ii]][3:7])
                    TIMEs = TIMEs.drop(labels=[ev_id_idx_copy[find_earliest_per_ii]])
                    ev_id_idx_copy.remove(ev_id_idx_copy[find_earliest_per_ii])
                #stations_arrival = [ii[3:7] for ii in event_catalogue[1]['station'][ev_id_idx]]
                Y = [None] * len(stations_arrival)
                Z = []
                for x in station_list:
                    if x in stations_arrival:
                        #print(x)
                        #print(stations_arrival.index(x),station_list.index(x))
                        Y[stations_arrival.index(x)] = station_list.index(x)
                        #print(Y)
                    else:
                        Z.append(station_list.index(x))
                Y = [ele for ele in Y if ele is not None]
                reordered_E = [E[ii] for ii in Y]
                rest_of_E = [E[ii] for ii in Z]
                reordered_E = reordered_E + rest_of_E
                reordered_N = [ii+1 for ii in reordered_E]
                reordered_Z = [ii+2 for ii in reordered_E]
                reordered_station_list_a = [station_list[ii] for ii in Y]
                reordered_station_list_b = [station_list[ii] for ii in Z]
                reordered_station_list = reordered_station_list_a + reordered_station_list_b
                trace_order_sh = 'event arrival time'
                
            earliest_arrival = UTCDateTime(ev_t1[ev_id_idx[earliest_arrival_idx]])
            LHS_t1 = earliest_arrival - arrival_buffer
            RHS_t2 = earliest_arrival + total_duration - arrival_buffer
            # Grab the stream between t1 - buffer_start to t2 + buffer_end per event
            stream_t1 = UTCDateTime(input_stream[0].stats.starttime)
            idx_starttime = int((LHS_t1 - stream_t1)*200)
            idx_endtime = int((RHS_t2 - stream_t1)*200)
    
    
            #########################################################################################
            #########################################################################################
            #########################################################################################
            # DRAWING TO PLOT:
                
            fig = plt.figure(figsize=(20, 9))
            
            # Draw the seismic signals, get the y-indices for station labelling
            yt = []
            icount = 0
            
            # Plotting ALL components stacked
            if component == 'ALL':
                for ii in range(len(reordered_E)): #42 stations and two line breaks
                    if '2C.%s.' %reordered_station_list[ii] in stations_that_det_event:
                        plt.plot(input_stream[reordered_E[ii]].data[idx_starttime:idx_endtime:sampling_interval] - vertical_scaling_range*(icount+1), color='k',linewidth=line_width)
                    else:
                        plt.plot(input_stream[reordered_E[ii]].data[idx_starttime:idx_endtime:sampling_interval] - vertical_scaling_range*(icount+1), color='w',alpha = 0.0,linewidth=line_width)                        
                    yt.append(0-vertical_scaling_range*(icount+1))
                    icount += 1
                for ii in range(len(reordered_N)): #42 stations and two line breaks
                    if '2C.%s.' %reordered_station_list[ii] in stations_that_det_event:
                        plt.plot(input_stream[reordered_N[ii]].data[idx_starttime:idx_endtime:sampling_interval] - vertical_scaling_range*(icount+1), color='k',linewidth=line_width)
                    else:
                        plt.plot(input_stream[reordered_N[ii]].data[idx_starttime:idx_endtime:sampling_interval] - vertical_scaling_range*(icount+1), color='w',alpha = 0.0,linewidth=line_width)                        
                    yt.append(0-vertical_scaling_range*(icount+1))
                    icount += 1
                for ii in range(len(reordered_Z)): #42 stations and two line breaks
                    if '2C.%s.' %reordered_station_list[ii] in stations_that_det_event:
                        plt.plot(input_stream[reordered_Z[ii]].data[idx_starttime:idx_endtime:sampling_interval] - vertical_scaling_range*(icount+1), color='k',linewidth=line_width)
                    else:
                        plt.plot(input_stream[reordered_Z[ii]].data[idx_starttime:idx_endtime:sampling_interval] - vertical_scaling_range*(icount+1), color='w',alpha = 0.0,linewidth=line_width)                        
                    yt.append(0-vertical_scaling_range*(icount+1))
                    icount += 1   
                    
                # Draw the shading per station
                icount = 0
                for ii in range(len(reordered_E)): #
                    if '2C.%s.' %reordered_station_list[ii] in stations_that_det_event:
                        # Get t1 of event start for each station
                        idx_ev = stations_that_det_event.index('2C.%s.' %reordered_station_list[ii])
                        idx_ev_starttime = int((UTCDateTime(list(ev_t1)[idx_ev]) - stream_t1)*200)/sampling_interval
                        idx_ev_endtime = max(int((UTCDateTime(list(ev_t1)[idx_ev])+list(ev_dur)[idx_ev] - stream_t1)*200),total_duration*200)/sampling_interval
                        plt.fill_between(np.arange(idx_ev_starttime-idx_starttime/sampling_interval,idx_ev_endtime+1-idx_starttime/sampling_interval,sampling_interval), yt[icount]-vertical_scaling_range, yt[icount]+vertical_scaling_range, facecolor=event_color, interpolate=True,alpha=event_alpha)       
                    icount += 1
                
                for ii in range(len(reordered_N)): #
                    if '2C.%s.' %reordered_station_list[ii] in stations_that_det_event:                
                    # Get t1 of event start for each station
                        idx_ev = stations_that_det_event.index('2C.%s.' %reordered_station_list[ii])
                        idx_ev_starttime = int((UTCDateTime(list(ev_t1)[idx_ev]) - stream_t1)*200)/sampling_interval
                        idx_ev_endtime = max(int((UTCDateTime(list(ev_t1)[idx_ev])+list(ev_dur)[idx_ev] - stream_t1)*200),total_duration*200)/sampling_interval
                        plt.fill_between(np.arange(idx_ev_starttime-idx_starttime/sampling_interval,idx_ev_endtime+1-idx_starttime/sampling_interval,sampling_interval), yt[icount]-vertical_scaling_range, yt[icount]+vertical_scaling_range, facecolor=event_color, interpolate=True,alpha=event_alpha)         
                    icount += 1
                    
                for ii in range(len(reordered_Z)): #
                    if '2C.%s.' %reordered_station_list[ii] in stations_that_det_event:
                    # Get t1 of event start for each station
                        idx_ev = stations_that_det_event.index('2C.%s.' %reordered_station_list[ii])
                        idx_ev_starttime = int((UTCDateTime(list(ev_t1)[idx_ev]) - stream_t1)*200)/sampling_interval
                        idx_ev_endtime = max(int((UTCDateTime(list(ev_t1)[idx_ev])+list(ev_dur)[idx_ev] - stream_t1)*200),total_duration*200)/sampling_interval
                        plt.fill_between(np.arange(idx_ev_starttime-idx_starttime/sampling_interval,idx_ev_endtime+1-idx_starttime/sampling_interval,sampling_interval), yt[icount]-vertical_scaling_range, yt[icount]+vertical_scaling_range, facecolor=event_color, interpolate=True,alpha=event_alpha)          
                    icount += 1  
                    
                plt.text(-1*arrival_buffer,yt[0],'HHE',color='black', 
                         bbox=dict(facecolor='white', edgecolor='black'),fontsize=12)
                plt.hlines(yt[0],-1*arrival_buffer,total_duration,colors='k')
                plt.text(-1*arrival_buffer,yt[14],'HHN',color='black', 
                         bbox=dict(facecolor='white', edgecolor='black'),fontsize=12)
                plt.hlines(yt[14],-1*arrival_buffer,total_duration,colors='k')
                plt.text(-1*arrival_buffer,yt[28],'HHZ',color='black', 
                         bbox=dict(facecolor='white', edgecolor='black'),fontsize=12)
                plt.hlines(yt[28],-1*arrival_buffer,total_duration,colors='k')
                plt.yticks(yt,reordered_station_list*3,fontsize=12)
            
            # IF PLOTTING JUST N COMPONENT
            elif component == 'N':    
                for ii in range(len(reordered_N)): #42 stations and two line breaks
                    if '2C.%s.' %reordered_station_list[ii] in stations_that_det_event:
                        plt.plot(input_stream[reordered_N[ii]].data[idx_starttime:idx_endtime:sampling_interval] - vertical_scaling_range*(icount+1), color='k',linewidth=line_width)
                    else:
                        plt.plot(input_stream[reordered_N[ii]].data[idx_starttime:idx_endtime:sampling_interval] - vertical_scaling_range*(icount+1), color='w',alpha = 0.0,linewidth=line_width)                        
                    yt.append(0-vertical_scaling_range*(icount+1))
                    icount += 1
                # Draw the shading per station
                icount = 0
                for ii in range(len(reordered_N)): #42 stations and two line breaks
                    if '2C.%s.' %reordered_station_list[ii] in stations_that_det_event:                
                    # Get t1 of event start for each station
                        idx_ev = stations_that_det_event.index('2C.%s.' %reordered_station_list[ii])
                        idx_ev_starttime = int((UTCDateTime(list(ev_t1)[idx_ev]) - stream_t1)*200)/sampling_interval
                        idx_ev_endtime = max(int((UTCDateTime(list(ev_t1)[idx_ev])+list(ev_dur)[idx_ev] - stream_t1)*200),total_duration*200)/sampling_interval
                        plt.fill_between(np.arange(idx_ev_starttime-idx_starttime/sampling_interval,idx_ev_endtime+1-idx_starttime/sampling_interval,sampling_interval), yt[icount]-vertical_scaling_range, yt[icount]+vertical_scaling_range, facecolor=event_color, interpolate=True,alpha=event_alpha)    
                    icount += 1        
                plt.text(-1*arrival_buffer,yt[0],'HHN',color='black', 
                         bbox=dict(facecolor='white', edgecolor='black'),fontsize=12)
                plt.yticks(yt,reordered_station_list,fontsize=12)
            
            # IF PLOTTING JUST E COMPONENT
            elif component == 'E':    
                for ii in range(len(reordered_E)): #42 stations and two line breaks
                    if '2C.%s.' %reordered_station_list[ii] in stations_that_det_event:
                        plt.plot(input_stream[reordered_E[ii]].data[idx_starttime:idx_endtime:sampling_interval] - vertical_scaling_range*(icount+1), color='k',linewidth=line_width)
                    else:
                        plt.plot(input_stream[reordered_E[ii]].data[idx_starttime:idx_endtime:sampling_interval] - vertical_scaling_range*(icount+1), color='w',alpha = 0.0,linewidth=line_width)                        
                    yt.append(0-vertical_scaling_range*(icount+1))
                    icount += 1
                # Draw the shading per station
                icount = 0
                for ii in range(len(reordered_E)): #42 stations and two line breaks
                    if '2C.%s.' %reordered_station_list[ii] in stations_that_det_event:                
                    # Get t1 of event start for each station
                        idx_ev = stations_that_det_event.index('2C.%s.' %reordered_station_list[ii])
                        idx_ev_starttime = int((UTCDateTime(list(ev_t1)[idx_ev]) - stream_t1)*200)/sampling_interval
                        idx_ev_endtime = max(int((UTCDateTime(list(ev_t1)[idx_ev])+list(ev_dur)[idx_ev] - stream_t1)*200),total_duration*200)/sampling_interval
                        plt.fill_between(np.arange(idx_ev_starttime-idx_starttime/sampling_interval,idx_ev_endtime+1-idx_starttime/sampling_interval,sampling_interval), yt[icount]-vertical_scaling_range, yt[icount]+vertical_scaling_range, facecolor=event_color, interpolate=True,alpha=event_alpha)     
                    icount += 1 
                plt.text(-1*arrival_buffer,yt[0],'HHE',color='black', 
                         bbox=dict(facecolor='white', edgecolor='black'),fontsize=12)
                plt.yticks(yt,reordered_station_list,fontsize=12)
            
            # IF PLOTTING JUST Z COMPONENT
            elif component == 'Z':    
                for ii in range(len(reordered_Z)): #42 stations and two line breaks
                    if '2C.%s.' %reordered_station_list[ii] in stations_that_det_event:
                        plt.plot(input_stream[reordered_Z[ii]].data[idx_starttime:idx_endtime:sampling_interval] - vertical_scaling_range*(icount+1), color='k',linewidth=line_width)
                    else:
                        plt.plot(input_stream[reordered_Z[ii]].data[idx_starttime:idx_endtime:sampling_interval] - vertical_scaling_range*(icount+1), color='w',alpha = 0.0,linewidth=line_width)                        
                    yt.append(0-vertical_scaling_range*(icount+1))
                    icount += 1
                # Draw the shading per station
                icount = 0
                for ii in range(len(reordered_Z)): #42 stations and two line breaks
                    if '2C.%s.' %reordered_station_list[ii] in stations_that_det_event:                
                    # Get t1 of event start for each station
                        idx_ev = stations_that_det_event.index('2C.%s.' %reordered_station_list[ii])
                        idx_ev_starttime = int((UTCDateTime(list(ev_t1)[idx_ev]) - stream_t1)*200)/sampling_interval
                        idx_ev_endtime = max(int((UTCDateTime(list(ev_t1)[idx_ev])+list(ev_dur)[idx_ev] - stream_t1)*200),total_duration*200)/sampling_interval
                        plt.fill_between(np.arange(idx_ev_starttime-idx_starttime/sampling_interval,idx_ev_endtime+1-idx_starttime/sampling_interval,sampling_interval), yt[icount]-vertical_scaling_range, yt[icount]+vertical_scaling_range, facecolor=event_color, interpolate=True,alpha=event_alpha)       
                    icount += 1     
                plt.text(-1*arrival_buffer,yt[0],'HHZ',color='black', 
                         bbox=dict(facecolor='white', edgecolor='black'),fontsize=12)
                plt.yticks(yt,reordered_station_list,fontsize=12)
    
            #########################################################################################
            # Get confidence marker if argument is True
            if confidence_marker == True:
                if len(stations_that_det_event) <= 5:
                    event_stations = []
                    for xx in range(len(stations_that_det_event)):
                        event_stations.append(stations_that_det_event[xx][3:7])
                    confidence_marker_code = confidence(event_stations,coincidence_sum,station_dict)
                else:
                    confidence_marker_code = 'High'
            else:
                confidence_marker_code = 'N/A'
            
            # _sh for shorthand
            highpass_sh = '%sHz' %(format(highpass,'.1E'))
            lowpass_sh = '%sHz' %(format(lowpass,'.1E'))
            
            ref_duration_in_HHMMSS = strftime("%H:%M:%S", gmtime(ref_dur))
            
            # Populate text box
            plt.text((total_duration*200 - 0.1*total_duration*200)/sampling_interval, 1*vertical_scaling_range, 't0 at %s \n Ref. event duration: %s \n %s confidence \n HP %s \n LP %s \n Order by: %s' %(UTCDateTime(earliest_arrival),\
                                                                                                 ref_duration_in_HHMMSS,confidence_marker_code,\
                                                                                                 highpass_sh,lowpass_sh,\
                                                                                                     trace_order_sh), color='black',
                                                                             bbox=dict(facecolor='white', alpha=0.2, edgecolor='black'),fontsize=10)
            
            plt.ylim(1*yt[len(yt)-1]-vertical_scaling_range,0)
            #plt.xticks(np.arange(0,(total_duration*200+200)/sampling_interval,(60*200)/sampling_interval),\
            #               ['-1','0','1','2','3','4','5','6','7','8','9'],fontsize=12)
            # Default xticks:
            plt.xticks(np.arange(0,(total_duration*200+200)/sampling_interval,(60*5*200)/sampling_interval),\
                           ['-5','0','5','10','15','20','25','30','35','40','45','50','55'],fontsize=12)
            plt.minorticks_on()
            ax = plt.gca()
            ax.xaxis.set_minor_locator(MultipleLocator(tick_lines_every_x_secs*200/sampling_interval))   # minor ticks every 10
            plt.xlabel('Time (minutes)',fontsize=12)             
            plt.xlim(0,(total_duration*200)/sampling_interval)
            plt.title('EVENT ID: %s [User label: %s]' %(ev_id,label),fontsize=22)
            #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
            plt.show()
            plt.pause(0.1)

            if interactive is False:
                pdf1.savefig(fig)
                plt.clf()
            elif interactive is True:
                plt.pause(0.1)
                print('User can now manipulate plot.')
                print('Reminder: User can manually save figure but \n figure will not save automatically in a .pdf file')
                zoom_ok = False
                print('\nZoom or pan to view, \npress spacebar when done (do NOT exit out of plot) \n')
                
                while not zoom_ok:
                    zoom_ok = plt.waitforbuttonpress()
                # Take an user input
                user_terminate = input("Only type Y+[Enter] in console when \n user is ready to move onto next plot: ")                    
                plt.close()
                        
        if interactive is False:    
            pdf1.close()
            plt.close('all')
        elif interactive is True and user_terminate != '':
            plt.close('all')
    
    return None

##################################################################################################
