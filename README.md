# RadarX
 
Python/NumPy port of Dr. Joey Bray's Radar-X 2019 Operation REDSHIELD - Transportable Ground Radar Engagement for EEE474. 

Written during the Summer of 2020 as I was learning how to use object-oriented programming and NumPy after losing a license for MATLAB. Now published (unfinished) to GitHub in June 2021.

radarProc.py introduces a Target class and an ADC class. The target class is essentially a struct to hold information about a target such as it's radial speed, elecation and SNR. Targets are created and then added to the ADC class. The ADC class keeps track of simulation time and locations of the targets (calculated upon each reading of the ADC). Reading the ADC will give an I/Q sample of the Pulse-Doppler Radar. Coherent processing intervals were not fully implemented.

As students we had to create a radar processor for a simulated pulse-Doppler radar. The processor had to disambiguate targets in range and speed (Doppler), use M-of-N binary integration, and, electronically steer a beam in azimuth. I have reimplemented some of this in the ADC_tests.ipynb notebook and included some of my group's original MATLAB scripts.


Radar Specifications:
| MODEL                        | CUSTOM AN/TPS                                                |
|------------------------------|--------------------------------------------------------------|
| TYPE                         | MEDIUM RANGE SURFACE/ARI PULSE DOPPLER                       |
| RADAR DIMENSIONS             | RANGE, AZIMUTH, ELEVATION (FIXED BEAMS), DOPPLER             |
| ELEVATION BEAMS              | 3 SELECTABLE (1 GROUND, 2 LOW ALTITUDE, 3 HIGH AIR)          |
| RF FREQUENCY                 | 6.0 GHz                                                      |
| PULSE COMPRESSION            | INACTIVE                                                     |
| PUSLE WIDTH                  | 1 MICROSECOND (FIXED)                                        |
| PEAK TRANSMIT RF POWER       | 10 KW                                                        |
| AZIMUTH FIELD OF VIEW        | 90 $\deg$ ($\pm 45\deg$ FROM BORESIGHT)                       |
| ARRAY ELEMENT SPACING        | $\lambda / 2$ IN AZIMUTH, 20 ELEMENTS IN AZIMUTH PER EL BEAM |
| AZIMUTH SCAN CONTROL         | SLECTABLE, PROVIDE PROGRESSIVE PHASE SHIFT $\Psi (\deg)$     |
| GAIN CONTROL                 | AUTOMATIC, AGC AND STC ENABLED                               |
| MAXIMUM DETECTION RANGE      | 30 KM                                                        |
| PULSE REPETITION FREQUENCIES | 500 TO 20 000 PULSES PER SECOND (SELECTABLE)                 |
| ANTENNA LOCATION             | GROUND                                                       |
| ADC OUTPUT                   | COHERENT DIGITAL I/Q                                         |
| SNR AT MAX RANGE             | 6 DB FOR A 1 SQUARE-METER TARGET, SINGLE PULSE               |
