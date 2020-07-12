# IntercomponentPhaseRelations
*Intercomponent Phase Relations* - framework for ICPR processing and analysis of statistics.

## Getting Started

Put your dataset into `In` folder. Put `comditions.m` function in datasets folder. It accepts file name and returns `conditions` structure contaning information about the current file. Specify dataset in `filePacketProcess.m` and start it. It will process each `*.mat` file and put ICPR data for each specified base frequency in each domain. Specify desirable ICPRs in `ICPRsuperPacket.m` and start it. See statistics at `Out/{your_dataset_name}`

## Common structure

* `/commonFunctions`		--------------- Service additional functions.
* `@signalObj`              --------------- Class for handling and processing real valued signal.
* `@quasiHarmonic`          --------------- Class for modelling harmonic and quasiharmonic signal.
* `@polyHarmonic`           --------------- Class for modelling polyharmonic signal.
* `@phaseProcessing`        --------------- Class for inter-component phase processing of polyharmonic signals.