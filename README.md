# Event detection for cryoseismology (Additional materials)
Latto, R., Turner, R. J., Reading, A. M., Winberry, J. P., “Event detection for cryoseismology.” The Cryosphere, in preparation for submission August 2022.

Here, users can find the data and code used for generating the event catalogues in Event detection for cryoseismology.
# Supplementary Materials: Chapter 3

Electronic Supplement Contents:
Supplementary files are publicly available on https://github.com/beccalatto/multistalta in sub folder
eventdetectionforcryoseismologyCh3.
Repository architecture, including folder and file names with descriptions, is listed below.

MultiSTALTAUseGuide.pdf: How-to guide to download the multi-STA/LTA function and use in
ObsPy
entrypoints.txt: Replacement text file for implementing multi-STA/LTA in Python
trigger.py: ObsPy script with all native trigger routines and including the multi-STA/LTA routine

## /prototype catalogues/detection catalogues/

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

## /prototype catalogues/msl labelled catalogues/

ReferenceEventCatalogueWhillansMultiSTALTAConfidenceAssignments.csv: Reference event
catalogue for the Whillans Ice Stream from the multi-STA/LTA algorithm with high and low confidence
assignments
ReferenceEventCatalogueWhillansMultiSTALTAKnownSeismicity.csv: Reference event cat-
alogue for the Whillans Ice Stream from the multi-STA/LTA algorithm with identified stick-slip and
teleseisms events
ReferenceEventCatalogueWhillansMultiSTALTATidesTable.txt: The .txt file contains all of

### 91


the reference catalogue columns as well as tidal height per event, tidal height derivative (increasing or
decreasing), absolute tidal height behavior (more or less positive, more or less negative), and inflection
type (minimum, maximum, or none). Tidal information results from the Circum-Antarctic Tidal Sim-
ulation for the location downstream from the grounding line of the Whillans Ice Stream (coordinate:
84 °20’20.3994”S -166°0’0”W) (Padman et al., 2002; Howard, 2019)

## /parameter evaluation/

Chapter3Section2.ipynb: Open access Python notebook used for the computational analysis and
compilation of figures included in Section 2 of Chapter 3. Notebook sections organized by...

```
□ 3.2.2 Simulation of test waveforms for algorithm development
```
```
□ 3.2.2.1 Representative event classes
□ 3.2.2.2 Simulated seismic signal
```
```
□ 3.2.3 Parameter search to optimize the application of the multi-STA/LTA
```
```
□ 3.2.3.1 Defining the fine-grid of parameters
□ 3.2.3.2 Assessing for probability of event detection
□ 3.2.3.3 Parameter search results
□ 3.2.3.5 Comparison of algorithms for synthetic data
```
## MyAnalystPlots/

MyAnalystPlots.py: Routine that plots events from the multistalta catalogue. The script enables
for event viewing in all three components (E,N,Z) and all stations that detect an event.
PDF products of the waveform views of the known Whillans Ice Stream seismicity, produced using the
MyAnalystPlots.py routine, made available separately due to file size at
https://cloudstor.aarnet.edu.au/plus/s/sLD2R6miP2wXoLY. Each product is divided into sets capped at
20 slides. The filenames are as follows:
MyAnalystPlotSTICK-SLIPPRATT14set0–6.pdf
MyAnalystPlotSTICK-SLIPPRATT14 ADDITIONAL.pdf
MyAnalystPlotTELESEISM Iset0–1.pdf
MyAnalystPlotTELESEISM IIset0–1.pdf
