# Multi-STA/LTA (Additional materials)

MultiSTALTA_UseGuide.pdf: How-to guide to download the multi-STA/LTA function and use in ObsPy

entrypoints.txt: Replacement text file for implementing multi-STA/LTA in Python

trigger.py: ObsPy script with all native trigger routines and including the multi-STA/LTA routine

<hr />

Repository architecture, including folder and file names with descriptions, is listed below.

# event_detection_for_cryoseismology/

 ## detection catalogues/

 ReferenceEventCatalogueWhillansMultiSTALTA.csv: Reference event catalogue for the Whillans
 Ice Stream from the multi-STA/LTA algorithm

 TraceEventCatalogueWhillansMultiSTALTA.csv : Trace event catalogue for the Whillans Ice
 Stream from the multi-STA/LTA algorithm

 ReferenceEventCatalogueWhillansRECmin.csv: Reference event catalogue for the Whillans Ice
 Stream from the recursive algorithm with RECmin parameters

 TraceEventCatalogueWhillansRECmin.csv: Trace event catalogue for the Whillans Ice Stream
 from the recursive algorithm with RECmin parameters

 ReferenceEventCatalogueWhillansRECmax.csv: Reference event catalogue for the Whillans Ice
 Stream from the recursive algorithm with RECmax parameters

 TraceEventCatalogueWhillansRECmax.csv: Trace event catalogue for the Whillans Ice Stream
 from the recursive algorithm with RECmax parameters

 ## labelled_catalogues/

 ReferenceEventCatalogueWhillansMultiSTALTAConfidenceAssignments.csv: Reference event
 catalogue for the Whillans Ice Stream from the multi-STA/LTA algorithm with high and low confidence
 assignments

 ReferenceEventCatalogueWhillansMultiSTALTAKnownSeismicity.csv: Reference event cat-
 alogue for the Whillans Ice Stream from the multi-STA/LTA algorithm with identified stick-slip and
 teleseisms events

 ReferenceEventCatalogueWhillansMultiSTALTATidesTable.txt: The .txt file contains all of the reference catalogue columns as well as tidal height per event, tidal height derivative (increasing or decreasing), absolute tidal height behavior (more or less positive, more or less negative), and inflection
 type (minimum, maximum, or none). Tidal information results from the Circum-Antarctic Tidal Sim-
 ulation for the location downstream from the grounding line of the Whillans Ice Stream (coordinate:
 84 °20’20.3994”S -166°0’0”W) (Padman et al., 2002; Howard, 2019)

 ## parameter_evaluation/

 EventDetectSectionS2.ipynb: Open access Python notebook used for the computational analysis and
 compilation of figures included in Section S2 of Event detection for cryoseismology. Notebook sections organized by...

   □ S2.2a: Simulation of test waveforms for algorithm development  <br />
       □ S2.2a.1: Representative event classes  <br />
       □ S2.2a.2: Simulated seismic signal  <br />
    
   □ S2.2b: Parameter search to optimize the application of the multi-STA/LTA  <br />
       □ S2.2b.1: Defining the fine-grid of parameters  <br />
       □ S2.2b.2: Assessing the success of event detection  <br />
       □ S2.2b.3: Parameter search results  <br />
       □ S2.2b.4: Parameter recommendations  <br />
       □ S2.2b.5: Comparison of algorithms for synthetic data  <br />

 ## MyAnalystPlots/

 MyAnalystPlots.py: Routine that plots events from the multistalta catalogue. The script enables
 for event viewing in all three components (E,N,Z) and all stations that detect an event.
 PDF products of the waveform views of the known Whillans Ice Stream seismicity, produced using the
 MyAnalystPlots.py routine, made available separately by request to author (due to file size). Each product is divided into sets capped at
 20 slides. The filenames are as follows:

 MyAnalystPlotSTICK-SLIPPRATT14set0–6.pdf

 MyAnalystPlotSTICK-SLIPPRATT14 ADDITIONAL.pdf

 MyAnalystPlotTELESEISM Iset0–1.pdf

 MyAnalystPlotTELESEISM IIset0–1.pdf

# unsupervised_learning_for_cryoseismology/

 TraceFeatureDatasetWhillans.csv: Trace feature dataset for the Whillans Ice Stream <br />
 ReferenceFeatureDatasetWhillans.csv: Reference feature dataset for the Whillans Ice Stream <br />
 ReferenceClusterLabelsWhillans k10.csv: Reference event catalogue with reference features and a 
 column for the cluster label from k-means++ applied to data from the Whillans Ice Stream, for k=10 

<hr />

# References:
Latto, R. (2022). Active Glacier Processes From Machine Learning Applied to Seismic Records [MSc Thesis]. University of Tasmania. School of Natural Sciences (Physics).

Latto, R., Turner, R. J., Reading, A. M., Winberry, J. P., “Event detection for cryoseismology.” The Cryosphere, in preparation for submission Feb 2024.

Latto, R., Turner, R. J., Reading, A. M., Cook, S., Kulessa, B., Winberry, J. P., “Unsupervised learning applied to cryoseismic signals: identification of glacier processes from the Whillans Ice Stream.” The Cryosphere, in preparation for submission Feb 2024.
