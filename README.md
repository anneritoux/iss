# How to run
The only file that you need to run to obtain and save the data is [bridge_process_template.m](https://github.com/jduffield65/iss/blob/master/bridge_process_template.m). The following will explain the changes to this file that need to made in order for it work with your data.

## File names
There are another of file/folder paths which need to be given:
* [```o.InputDirectory```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L7): This is the path of the folder that contains the raw data, of type specified by [```o.RawFileExtension```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L20), e.g. ```.nd2```. An example of what this folder typically looks like is given below:

<p float="left">
<img src="DebugImages/README/InputDirectory.png" width = "450"> 
</p>

* [```o.FileBase```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L10-L18): These are the names of the files within ```o.InputDirectory``` (minus the extension). For the above example, we would set
<pre>
o.FileBase{1} = 'Exp1_r0';
o.FileBase{2} = 'Exp1_r1';
&#8942
o.FileBase{7} = 'Exp1_r6';
o.FileBase{8} = 'Exp1_anchor';
</pre>
You need to make sure that ```o.FileBase{```[```o.ReferenceRound```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L24)```}``` is set to the anchor round. Also, the other rounds must be in the correct imaging order.

* [```o.TileDirectory```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L21): This is the path for the folder that you would like the filtered images for each tile, round and colour channel to be saved to. The file named as  ```o.FileBase{r}_tT.tif``` contains all the colour channel images for round r, tile T.

* [```o.OutputDirectory```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L26): This is the path for the folder that you would like the iss objects after each step of the pipeline to be saved.

* [```o.CodeFile```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L105): This is the path for the file containing the code for each gene. The file should be a text file containing two columns, the first being the gene name. The second is the code of length [```o.nRounds```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L5), containing numbers in the range from 0 to [```o.nBP```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L74)-1. An example codebook for ```o.nRounds=o.nBP=7``` is shown below:

<p float="left">
<img src="DebugImages/README/CodeBook.png" width = "250"> 
</p>

## Stitching and registration parameters
These are parameters that slightly affect how the stitching of tiles and registration between rounds and colour channels works. The default values in [bridge_process_template.m](https://github.com/jduffield65/iss/blob/master/bridge_process_template.m) should work most of the time but there are some cases when they may need altering.

* [```o.InitialShiftChannel```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L82): This is the channel used to register all rounds to the anchor, so ensure it is one of the best colour channels. 

* [```o.RegSearch```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L61-L64): The algorithm for stitching together tiles only looks at shifts in the range specified by ```o.RegSearch```. The values that may need changing here are ```o.RegSearch.South.Y``` and ```o.RegSearch.East.X```. The default values of ```-1900:o.RegStep(2):-1700``` are heavily dependent on [```o.TileSz```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L44). I.e. we only consider overlaps between 148 - 348 pixels which covered the expected value of about 10% (205 pixels). If the expected overlap or [```o.TileSz```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L44) is different to this though, these values will need changing. 

* [```o.FindSpotsSearch```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L88-L89): The algorithm for finding shifts between the anchor round and each imaging round only looks at shifts in the range specified by ```o.FindSpotsSearch```. The default values assume a relatively small shift between of absolute value less than 100 pixels. However, if you know one round had a particularly large shift, you will need to change this range. You can also specify a different range for each round through:
```matlab
o.FindSpotsSearch = cell(o.nRounds,1);
for r = o.UseRounds
    o.FindSpotsSearch{r}.Y = MinYShift:o.FindSpotsStep(1):MaxYShift;
    o.FindSpotsSearch{r}.X = MinXShift:o.FindSpotsStep(1):MaxXShift;
end
```

## Visualising results
The pipeline runs two different algorithms for assigning spots to genes. The following will give instructions explaining how to view the distribution of assigned genes in both cases.

### Dot product method
To visualise the results, load in the final saved iss object which should be named [```oCallSpots.mat```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L111). Then, load in the  background dapi image and run [```o.plot```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L116-L117). This will show you the gene assignments (saved as [```o.SpotCodeNo```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L519)) given by the file [```o.call_spots```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L109) which is achieved by taking the dot product of the [normalised spot](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L552) and [gene codes](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L551). Only the results where this dot product (saved as [```o.SpotScore```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L524)) is above [```o.CombiQualThresh```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L115) are shown. An example plot is given below with ```o.CombiQualThresh = 0.7```.

<p float="left">
<img src="DebugImages/README/CombiThresh0.7.png" width = "450"> 
</p>

To change the value of the threshold, simply set ```o.CombiQualThresh = NewValue``` and then run ```iss_change_plot(o)```. The above data with ```o.CombiQualThresh = 0.9``` is shown below:

<p float="left">
<img src="DebugImages/README/CombiThresh0.9.png" width = "450"> 
</p>

With [```o.CallSpotsCodeNorm = 'WholeCode'```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L108), each spot and gene code has L2 norm of 1 so the maximum value of the dot product hence ```o.CombiQualThresh``` is 1 (Recommend ```o.CombiQualThresh ~ 0.7```). However, with ```o.CallSpotsCodeNorm = 'Round'```, each round in each code has L2 norm of 1 so each code has L2 norm of [```o.nRounds```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L5). So in this case, the max value of the dot product hence ```o.CombiQualThresh``` would be ```o.nRounds``` (Recommend ```o.CombiQualThresh ~ 4``` for ```o.nRounds=7```. Note that if you change ```o.CallSpotsCodeNorm```, then you need to run [```o.call_spots```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L109) again. The justification for setting ```o.CallSpotsCodeNorm = 'Round'``` is that with no bleed through, we expect each spot to appear in one colour channel in each round so we want to give each round equal weighting, but with ```o.CallSpotsCodeNorm = 'WholeCode'```, a particularly intense round would dominate the others.

### Probability method
To view the results in this case, start off as before by loading in the Dapi image and running [```o.plot```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L116-L117). Then to see the probability gene assignments, run [```iss_change_plot(o,'Prob')```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L122). These are the gene assignments (saved as [```o.pSpotCodeNo```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L652)) given by [```o.call_spots_prob```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L110)). 

This method works by finding the probability, in each round and channel, that the result of removing a scaled version of the gene code from the spot code can be explained by the background distribution in that round and channel. Then [total log probability](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L642) is found by summing up all ```o.nRounds*o.nBP``` of the log of these probabilities. The score used for this method (saved as [```o.pSpotScore```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L645)) is then the largest log probability minus the second largest i.e. the probability the spot is the gene given by ```o.pSpotCodeNo``` minus the probability that the spot is the next most likely gene. 

Spots for which ```o.pSpotScore>o.pScoreThresh``` or ```o.pSpotIntensity>o.pIntensityThresh``` are then shown. For spot s, [```o.pSpotIntensity(s)```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L639) is the mean spot intensity of the ```o.nRounds``` spot intensities specified by [```o.CharCodes(o.pSpotNo(s))```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L536) minus the mean of all of the other ```o.nRounds*o.nBP-o.nRounds``` spot intensities in the code. This is used as the probability method tends to penalise the score if the intensity of the spot is higher than expected by the gene code. This feature is not desirable, but ```o.pSpotIntensity``` would benefit from such anomalously large intensities but only if they are in the rounds/channels predicted by the gene it was assigned to. The plot for the same data set shown for the dot productmethod is given below, with ```o.pScoreThresh=10``` and ```o.pIntensityThresh=100``` (I would recommend to use values near these).

<p float="left">
<img src="DebugImages/README/Score10Intensity100.png" width = "450"> 
</p>

As with the other method, to change the value of the threshold, simply set ```o.pScoreThresh = NewValue``` and  ```o.pIntensityThresh = NewValue```and then run ```iss_change_plot(o)```. The above data with ```o.pScoreThresh = 20``` and ```o.pIntensityThresh = 1000``` is shown below:

<p float="left">
<img src="DebugImages/README/Score20Intensity1000.png" width = "450"> 
</p>

### Which method to use?
The dot product method involves relative normalisation between rounds and colour channels to make them more equal and thus have a more equal contribution to the dot product. However, this sometimes causes the worse colour channels (usually one and three) to be boosted too much causing false assignments. An example of this is given below (codes and spots are normalised so have L2 norm of 1).

Here, because round 7, channel 6  is particularly low intensity, when it is normalised it gets boosted resulting in this square dominating the whole code. Then to have a high dot product, this spot must match to a gene which is also high in round 7, channel 6 even though it doesn't match any other squares.


<p float="left">
<img src="DebugImages/README/BadDotProduct.png" width = "450"> 
</p>

The probability method does not involve any such normalisation so is probably the better method to use. Also, with the recommended threshold values (first plot of each method), much more spots overcome the thresholding than the dot product method (45% more in this example).

### Visualising individual spots
To view the dot product assignment of a particular gene, with the plot open, run [```iss_view_codes(o,234321,Norm)```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L119). This will cause a crosshair to appear on the plot, then just click on the spot of interest as shown below.

<p float="left">
<img src="DebugImages/README/Crosshair.png" width = "450"> 
</p>

234321 is just the figure number of the plot (should always be the same). Norm controls the normalisation applied to the spot and gene codes. You can set Norm equal to 1,2 or 3 to highlight certain features:
* Norm = 1: This gives the raw values, for the spot in the previous section, this would be:

<p float="left">
<img src="DebugImages/README/Norm1.png" width = "450"> 
</p>

* Norm = 2: This normalises in the same way used in ```call_spots``` i.e. it normalises by the percentile given by [```o.SpotNormPrctile```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L267) in each colour channel and round and then depending on ```o.CallSpotsCodeNorm```, it normalises the resultant so the either the whole code has L2 norm of 1 or each round does. The plot in the previous section used Norm = 2.

* Norm = 3: This normalises each colour channel by the percentile given by ```o.SpotNormPrctile```  across all rounds. Using this, the spot we are considering would appear like this:

<p float="left">
<img src="DebugImages/README/Norm3.png" width = "450"> 
</p>

Also, if you know the index (```SpotNo```) of the spot you are interested in but don't want to find it in the plot, you can just run ```iss_view_codes(o,234321,Norm,SpotNo)```.

The equivalent plot for the probability method is [```iss_view_prob(o,234321,Norm)```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L123) which is called in exactly the same way. The plot is the same except the gene shown is that given by ```o.pSpotCodeNo``` not ```o.SpotCodeNo```. Also, there is an extra plot which shows you how the overall probability is made up from the individual probabilities of each round and channel. This plot for the example spot is given below:

<p float="left">
<img src="DebugImages/README/view_prob.png" width = "450"> 
</p>

The plot also gives the spot score, log probability, score deviation and intensity values for its assignment. LogProb is the total probability i.e. sum of all squares in bottom plot. Score deviation is the standard deviation of the log probabilities for assigning the spot to every gene in the codebook. This is included as if ```o.pSpotScoreDev+o.pSpotScore<o.pDevThresh```, then the assignment is rejected. This is to get rid of spots that have a similar probability when matched to every gene. Also, the values of these variables are coloured - green means that value caused the gene to be automatically accepted, red means that thay value caused the gene to fail (green supersedes red most of the time). If all variables are black, then the match is also accepted.

### Understanding the probability method
If you left/right click on a particular square in the bottom plot of the ```iss_view_prob``` window, the plot on the left/right below will appear further explaining the probability for that particular round and channel. These particular plots are for round 7, channel 2 for the above example spot.

<p float="left">
<img src="DebugImages/README/NonCodeGoodPlt1.png" width = "400"> 
<img src="DebugImages/README/NonCodeGoodPlt2.png" width = "400">
</p>
