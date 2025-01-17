# How to run
The only file that you need to run to obtain and save the data is [bridge_process_template.m](https://github.com/jduffield65/iss/blob/master/bridge_process_template.m). The following will explain the changes to this file that need to made in order for it work with your data.

## Parameters that should be checked before each run
There are a few parameters that need double checking or adjusting before each run:
* [```o.AnchorRound```](https://github.com/jduffield65/iss/blob/b537681136244984efc1182d23f244a7e3dc9caf/bridge_process_template.m#L4): This is the index of the anchor round. Should be the first round after imaging rounds. The anchor round is the round that contains the Dapi image and anchor channel image, which contains spots arising from all genes.
* [```o.AnchorChannel```](https://github.com/jduffield65/iss/blob/b537681136244984efc1182d23f244a7e3dc9caf/bridge_process_template.m#L5): This is the channel (starting from 1) in the anchor round (given by [```o.AnchorRound```](https://github.com/jduffield65/iss/blob/b537681136244984efc1182d23f244a7e3dc9caf/bridge_process_template.m#L4)) which contains the spots. The final genes (in probability and dot product methods) stem from these spots and it is also used for registration.
* [```o.DapiChannel```](https://github.com/jduffield65/iss/blob/b537681136244984efc1182d23f244a7e3dc9caf/bridge_process_template.m#L6): This is the channel in the anchor round that contains the Dapi images.
* [```o.InitialShiftChannel```](https://github.com/jduffield65/iss/blob/b537681136244984efc1182d23f244a7e3dc9caf/bridge_process_template.m#L7): This is the channel used to register all rounds to the anchor, so ensure it is one of the best colour channels. 
* [```o.ReferenceRound```](https://github.com/jduffield65/iss/blob/b537681136244984efc1182d23f244a7e3dc9caf/bridge_process_template.m#L8): This is the index of the reference round. The reference round, is the round that the global coordinate system is built upon and all other rounds are registered to, it can be equal to or different to [```o.AnchorRound```](https://github.com/jduffield65/iss/blob/b537681136244984efc1182d23f244a7e3dc9caf/bridge_process_template.m#L4). If it is not equal to [```o.AnchorRound```](https://github.com/jduffield65/iss/blob/b537681136244984efc1182d23f244a7e3dc9caf/bridge_process_template.m#L4), then the anchor round is not used for anything.
* [```o.ReferenceChannel```](https://github.com/jduffield65/iss/blob/b537681136244984efc1182d23f244a7e3dc9caf/bridge_process_template.m#L9): This is the channel (starting from 1) within the reference round (given by [```o.ReferenceRound```](https://github.com/jduffield65/iss/blob/b537681136244984efc1182d23f244a7e3dc9caf/bridge_process_template.m#L8)) which the global coordinate system is built upon. Ensure it is one of the best colour channels.
* [```o.RawFileExtension```](https://github.com/jduffield65/iss/blob/b537681136244984efc1182d23f244a7e3dc9caf/bridge_process_template.m#L10): This is the format of the raw data in [```o.InputDirectory```](https://github.com/jduffield65/iss/blob/f739d4c2e38c66ff82e7fd7a9f02b0fe73125353/bridge_process_template.m#L13).
* [```o.LogToFile```](https://github.com/jduffield65/iss/blob/5ba65f7c264798e66a91417be886707760062958/bridge_process_template.m#L11): Set this to 1 if you want the contents of the command window output to a .txt file. If you don't want the file, set to 0. By default, it is 1.

## File names
There are a number of file/folder paths which need to be given:
* [```o.InputDirectory```](https://github.com/jduffield65/iss/blob/f739d4c2e38c66ff82e7fd7a9f02b0fe73125353/bridge_process_template.m#L13): This is the path of the folder that contains the raw data, of type specified by [```o.RawFileExtension```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L20), e.g. ```.nd2```. An example of what this folder typically looks like is given below:

<p float="left">
<img src="DebugImages/README/InputDirectory.png" width = "450"> 
</p>

* [```o.FileBase```](https://github.com/jduffield65/iss/blob/f739d4c2e38c66ff82e7fd7a9f02b0fe73125353/bridge_process_template.m#L15-L24): These are the names of the files within ```o.InputDirectory``` (minus the extension). For the above example, we would set
<pre>
o.FileBase{1} = 'Exp1_r0';
o.FileBase{2} = 'Exp1_r1';
&#8942
o.FileBase{7} = 'Exp1_r6';
o.FileBase{8} = 'Exp1_anchor';
</pre>
You need to make sure that ```o.FileBase{```[```o.ReferenceRound```](https://github.com/jduffield65/iss/blob/f739d4c2e38c66ff82e7fd7a9f02b0fe73125353/bridge_process_template.m#L8)```}``` is set to the anchor round. Also, the other rounds must be in the correct imaging order.

* [```o.TileDirectory```](https://github.com/jduffield65/iss/blob/f739d4c2e38c66ff82e7fd7a9f02b0fe73125353/bridge_process_template.m#L26): This is the path for the folder that you would like the filtered images for each tile, round and colour channel to be saved to. The file named as  ```o.FileBase{r}_tT.tif``` contains all the colour channel images for round r, tile T.

* [```o.OutputDirectory```](https://github.com/jduffield65/iss/blob/f739d4c2e38c66ff82e7fd7a9f02b0fe73125353/bridge_process_template.m#L27): This is the path for the folder that you would like the iss objects after each step of the pipeline to be saved.

* [```o.CodeFile```](https://github.com/jduffield65/iss/blob/f739d4c2e38c66ff82e7fd7a9f02b0fe73125353/bridge_process_template.m#L30): This is the path for the file containing the code for each gene. The file should be a text file containing two columns, the first being the gene name. The second is the code specifying which colour channel that gene should appear in each round. Thus it is of length [```o.nRounds```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L5), containing numbers in the range from 0 to [```o.nBP```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L74)-1. An example codebook for ```o.nRounds=o.nBP=7``` is shown below:

<p float="left">
<img src="DebugImages/README/CodeBook.png" width = "250"> 
</p>

## Running with a subset of tiles
Running the pipeline with the whole set of tiles can be quite time consuming so it is sometimes desirable to run the pipeline on a subset of tiles. You can do this via the [```o.EmptyTiles```](https://github.com/jduffield65/iss/blob/9b863b1ff3589794334479cad0f31ce3db3698e3/%40iss/iss.m#L455-L456) variable. For example, to only use tiles 1 and 2, you can add the [following lines](https://github.com/jduffield65/iss/blob/9b863b1ff3589794334479cad0f31ce3db3698e3/bridge_process_template.m#L62-L65) to [bridge_process_template.m](https://github.com/jduffield65/iss/blob/PixelBased/bridge_process_template.m):

```matlab
o.EmptyTiles(:) = 1;
UseTiles = [1,2];
o.EmptyTiles(UseTiles) = 0;    %If after extract_and_filter
% o.EmptyTiles = UseTiles      %If before extract_and_filter
```

All tiles ```t```, such that ```o.EmptyTiles(t) = 1``` will be skipped. If you add it before the [registration step](https://github.com/jduffield65/iss/blob/9b863b1ff3589794334479cad0f31ce3db3698e3/bridge_process_template.m#L99), then if more than one tile is specified, they should each have at least one neighbour. An example showing three valid entries and one incorrect entry of ```o.EmptyTiles```, for a dataset consisting of 6 tiles is shown below.

:heavy_check_mark: | :heavy_check_mark: |  :heavy_check_mark: | :x:
:------------ | :-------------| :-------------| :-------------
| ```1 1```       | ```0 1```     |   ```0 0```   |  ```0 1```      
| ```0 0```       | ```0 1```     |   ```0 1```   |  ```1 0```      
| ```1 1```       | ```1 1```     |   ```1 1```   |  ```1 1```     
|Tiles 2 and 5 | Tiles 1 and 2 |  Tiles 1,2 and 4 | Tiles 1 and 5

The different line needed to specify tiles before extract_and_filter is because the algorithm hasn't yet worked out the location of the tiles.

Running the full pipeline (post extract_and_filter) should take on the order of half an hour, if only one tile is selected.

## Stitching and registration parameters
These are parameters that slightly affect how the stitching of tiles and registration between rounds and colour channels works. The default values in [bridge_process_template.m](https://github.com/jduffield65/iss/blob/master/bridge_process_template.m) should work most of the time but there are some cases when they may need altering.

* [```o.RegSearch```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L61-L64): The algorithm for stitching together tiles only looks at shifts in the range specified by ```o.RegSearch```. The values that may need changing here are ```o.RegSearch.South.Y``` and ```o.RegSearch.East.X```. The default values of ```-1900:o.RegStep(2):-1700``` are heavily dependent on [```o.TileSz```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L44). I.e. we only consider overlaps between 148 - 348 pixels which cover the expected value of about 10% (205 pixels). If the expected overlap or [```o.TileSz```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L44) is different to this though, these values will need changing. 

* [```o.FindSpotsSearch```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L88-L89): The algorithm for finding shifts between the anchor round and each imaging round only looks at shifts in the range specified by ```o.FindSpotsSearch```. The default values assume a relatively small shift of absolute value less than 100 pixels. However, if you know one round had a particularly large shift, you will need to change this range. You can also specify a different range for each round through:
```matlab
o.FindSpotsSearch = cell(o.nRounds,1);
for r = o.UseRounds
    o.FindSpotsSearch{r}.Y = MinYShift:o.FindSpotsStep(1):MaxYShift;
    o.FindSpotsSearch{r}.X = MinXShift:o.FindSpotsStep(1):MaxXShift;
end
```

## Loading old data and changing classes
The master branch of the pipeline has been changed so there are different classes relating to the different gene calling algorithms described in the next section. As a result, iss objects saved using an older version cannot be loaded into this branch directly.

To use this branch with old data, you can do the fllowing:
* Add both the new master branch and the old branch (the version of the iss object that the old data was run with) to the path.
* Load in the old data
* Run [```o = convert_iss_data(o);```](https://github.com/jduffield65/iss/blob/1dab53fdecf5312aeae33034a50061bcb7fb44e9/convert_iss_data.m)

This will convert it to an object of the class [iss_PixelBased](https://github.com/jduffield65/iss/blob/master/%40iss_PixelBased/iss_PixelBased.m). If you then want to run the OMP algorithm on this data, you will need to change the class to one of the OMP classes. To change class, you run [```o = switch_class(o, NewClass);```](https://github.com/jduffield65/iss/blob/1dab53fdecf5312aeae33034a50061bcb7fb44e9/@iss_Base/switch_class.m) where NewClass can be ```iss_Base, iss_GroundTruth, iss_OMP, iss_OMP, iss_OMP_ConstantBackground_WeightDotProduct, iss_PixelBased``` or ```iss_Spatial```. If the class change is such that data will be lost, an error will stop it from happening. If you still want to proceed, you can run ```o = switch_class(o, NewClass, true);```.

So, to run the OMP algorithm with old data you run:

```matlab
o = convert_iss_data(o);
o = switch_class(o, iss_OMP);
o = o.call_spots_omp;
```


## Visualising results
The pipeline can run three different algorithms for assigning genes to spots. The following will give instructions explaining how to view the distribution of assigned genes in each case.

### Dot product method
To visualise the results, load in the final saved iss object which should be named [```oCallSpots.mat```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L111). Then, load in the  background dapi image and run [```o.plot```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L116-L117). This will show you the gene assignments (saved as [```o.SpotCodeNo```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L519)) given by the file [```o.call_spots```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L109) which is achieved by taking the dot product of the [normalised spot](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L552) and [gene codes](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L551). Only the results where this dot product (saved as [```o.SpotScore```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L524)) is above [```o.CombiQualThresh```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L115) will be shown. An example plot is given below with ```o.CombiQualThresh = 0.7```.

<p float="left">
<img src="DebugImages/README/CombiThresh0.7.png" width = "650"> 
</p>

To change the value of the threshold, simply set ```o.CombiQualThresh = NewValue``` and then run ```iss_change_plot(o)```. The above data with ```o.CombiQualThresh = 0.9``` is shown below:

<p float="left">
<img src="DebugImages/README/CombiThresh0.9.png" width = "650"> 
</p>

With [```o.CallSpotsCodeNorm = 'WholeCode'```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L108), each spot and gene code has L2 norm of 1 so the maximum value of the dot product hence ```o.CombiQualThresh``` is 1 (Recommend ```o.CombiQualThresh ~ 0.7```). However, with ```o.CallSpotsCodeNorm = 'Round'```, each round in each code has L2 norm of 1 so each code has L2 norm of [```o.nRounds```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L5). So in this case, the max value of the dot product hence ```o.CombiQualThresh``` would be ```o.nRounds``` (Recommend ```o.CombiQualThresh ~ 4``` for ```o.nRounds=7```). Note that if you change ```o.CallSpotsCodeNorm```, then you need to run [```o.call_spots```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L109) again. The justification for setting ```o.CallSpotsCodeNorm = 'Round'``` is that with no bleed through, we expect each spot to appear in one colour channel in each round so we want to give each round equal weighting, but with ```o.CallSpotsCodeNorm = 'WholeCode'```, a particularly intense round would dominate the others.

### Probability method
To view the results in this case, start off as before by loading in the Dapi image and running [```o.plot```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L116-L117). Then to see the probability gene assignments, run [```iss_change_plot(o,'Prob')```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L122). These are the gene assignments (saved as [```o.pSpotCodeNo```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L652)) given by [```o.call_spots_prob```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L110)). 

This method works by finding the probability, in each round and channel, that the result of removing a scaled version of the gene code from the spot code can be explained by the background distribution in that round and channel. Then the [total log probability](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L642) is found by summing up all ```o.nRounds*o.nBP``` of the log of these probabilities. This is explained further in the section [Understanding the probability method](#understanding-the-probability-method). The score used for this method (saved as [```o.pSpotScore```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L645)) is then the largest log probability minus the second largest i.e. the probability the spot is the gene given by ```o.pSpotCodeNo``` minus the probability that the spot is the next most likely gene. 

Spots for which ```o.pSpotScore>o.pScoreThresh``` or ```o.pSpotIntensity>o.pIntensityThresh``` are then shown. For spot s, [```o.pSpotIntensity(s)```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L639) is the mean spot intensity of the ```o.nRounds``` spot intensities specified by [```o.CharCodes(o.pSpotNo(s))```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L536) minus the mean of all of the other ```o.nRounds*o.nBP-o.nRounds``` spot intensities in the code. This is used as the probability method tends to penalise the score if the intensity of the spot is higher than expected by the gene code. This feature is not desirable, but ```o.pSpotIntensity``` would benefit from such anomalously large intensities but only if they are in the rounds/channels predicted by the gene it was assigned to. The plot for the same data set shown for the dot product method is given below, with ```o.pScoreThresh=10``` and ```o.pIntensityThresh=100``` (I would recommend to use values near these).

<p float="left">
<img src="DebugImages/README/Score10Intensity100.png" width = "650"> 
</p>

As with the other method, to change the value of the threshold, simply set ```o.pScoreThresh = NewValue``` and  ```o.pIntensityThresh = NewValue```and then run ```iss_change_plot(o)```. The above data with ```o.pScoreThresh = 20``` and ```o.pIntensityThresh = 1000``` is shown below:

<p float="left">
<img src="DebugImages/README/Score20Intensity1000.png" width = "650"> 
</p>

### Pixel based results
To view the [results from the pixel based method](https://github.com/jduffield65/iss/blob/849350e6f0a4742d8fd6a3e083b0ffcd81914e31/%40iss/iss.m#L743-L771), run ```iss_change_plot(o,'Pixel')```. In this case, the spots are detected differentely, but the gene assignments are still carried out using the probability method. The thresholds to use are thus: ```o.pScoreThresh``` and ```o.pIntensityThresh```. The results saved have analagous names and meanings as with the probability method, except the prefix is ```px``` instead of ```p``` e.g. ```o.pxSpotScore```.

Also, the pixel based method allows for the possibility of multiple genes assigned to the same pixel. To view these overlapping genes, you can set ```o.pScoreThresh2``` to a value below 0. It has a default value of 0 meaning only genes that are the best match at each pixel can be shown. If you set it to ```o.pScoreThresh2 = -0.001;```, then it allows for spots for which ```o.pxSpotScore = 0``` i.e. the second best match at that pixel.

### OMP (Orthogonal Matching Pursuit) method
Start off as before, running [```o.plot```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L116-L117). Then run ```iss_change_plot(o,'OMP')```. These are the gene assignments (saved as [```o.ompSpotCodeNo```](https://github.com/jduffield65/iss/blob/4e0d03d53ad006c92073db3f20b9f6fd21557f0a/%40iss_OMP/iss_OMP.m#L150-L151)) given by [```o.call_spots_omp```](https://github.com/jduffield65/iss/blob/4e0d03d53ad006c92073db3f20b9f6fd21557f0a/bridge_process_template.m#L164).

This method works by performing an orthogonal matching pursuit algorithm on each pixel separately. First some [background vectors](https://github.com/jduffield65/iss/blob/4e0d03d53ad006c92073db3f20b9f6fd21557f0a/@iss_OMP/get_background_codes.m) [are fit](https://github.com/jduffield65/iss/blob/4e0d03d53ad006c92073db3f20b9f6fd21557f0a/omp_free_background.m#L49-L51) to explain the non-gene variation in the pixel. These are usually just 7 codes, each being a strip in a colour channel, as shown below.

<p float="left">
<img src="DebugImages/README/ompBackgroundVectors.png" width = "450"> 
</p>

Next, a [gene is selected](https://github.com/jduffield65/iss/blob/4e0d03d53ad006c92073db3f20b9f6fd21557f0a/omp_free_background.m#L61-L65) that best explains the [residual](https://github.com/jduffield65/iss/blob/4e0d03d53ad006c92073db3f20b9f6fd21557f0a/omp_free_background.m#L69) (pixel signal once background removed). Next, the [coefficient of this gene is found](https://github.com/jduffield65/iss/blob/4e0d03d53ad006c92073db3f20b9f6fd21557f0a/omp_free_background.m#L68) (i.e. how intense this gene is in this pixel). If the difference in the L2 norm of the residual before and after fitting the gene is [less than a threshold](https://github.com/jduffield65/iss/blob/4e0d03d53ad006c92073db3f20b9f6fd21557f0a/omp_free_background.m#L71-L78), the gene is rejected and we say this pixel contains no genes. If it [exceeds the threshold](https://github.com/jduffield65/iss/blob/4e0d03d53ad006c92073db3f20b9f6fd21557f0a/omp_free_background.m#L79-L83), the gene is accepted and a new iteration starts, fitting the best gene that can explain the new residual. This process continues until the residual difference falls below the threshold. At each step of the iteration, the coefficients of the previously added genes are also updated to account for the new gene. 

The OMP plot shows spots which pass quite a complicated thresholding process as described in the document [Thresholding.md](https://github.com/jduffield65/iss/blob/a8c4d104274775dae67c2bace447476f40f0f355/@iss_OMP_ConstantBackground_WeightDotProduct/Thresholding.md). Thresholding is done on three variables: ```o.ompNeighbNonZeros```, ```o.ompSpotIntensity``` and ```o.ompSpotScore```. There are three thresholds for each of these parameters. Using the default thresholds, the plot is the one below.

<p float="left">
<img src="DebugImages/README/ompResultsDefault.jpg" width = "650"> 
</p>

The threshold to alter to see the most obvious change is ```o.ompIntensityThresh2```. Increasing this from 0.001 to 0.1 and then running ```iss_change_plot(o)``` gives the following plot:

<p float="left">
<img src="DebugImages/README/ompResultsIncreaseIntensityThresh.jpg" width = "650"> 
</p>


### Which method to use?
The dot product method involves relative normalisation between rounds and colour channels to make them more equal and thus have a more equal contribution to the dot product. However, this sometimes causes the worse colour channels (usually one and three) to be boosted too much, causing false assignments. An example of this is given below (codes and spots are normalised so have L2 norm of 1).

<p float="left">
<img src="DebugImages/README/BadDotProduct.png" width = "450"> 
</p>

Here, because round 7, channel 6  is particularly low intensity, when it is normalised it gets boosted resulting in this square dominating the whole code. Then to have a high dot product, this spot must match to a gene which is also high in round 7, channel 6 even though it doesn't match any other squares.

The probability method does not involve any such normalisation so is probably the better method to use. The main advantages of the probability method is that it uses the actual measured background distribution and it allows for variation in intensity between rounds i.e. just because spot is intense in one round doesn't mean it is intense in all other rounds. However, it has a problem that it doesn't really allow for overlapping spots. 

The OMP method does allow for overlapping spots by fitting multiple genes to each pixel. It doesn't use the measured background distribution but the use of the background vectors helps get around this issue in practise. Neither does it allow for variation in intensity between rounds for a particular gene as just a single coefficient is found for each gene. To get around this, we [find some initial spots for each gene](https://github.com/jduffield65/iss/blob/a8c4d104274775dae67c2bace447476f40f0f355/@iss_OMP/call_spots_omp_initial.m) and use these to find the mean bled code for each gene. Using these mean codes, we can determine how [intense each gene is expected to be in each round compared to the bleed matrix prediction](https://github.com/jduffield65/iss/blob/a8c4d104274775dae67c2bace447476f40f0f355/@iss_OMP/get_gene_efficiencies.m). This gives us an updated codebook which accounts for round to round and gene to gene intensity variation. For example, the plot below shows that the OMP method accounts for the fact that Aldoc is especially weak in round 5:

<p float="left">
<img src="DebugImages/README/ompGeneEfficiencyExample.png" width = "450"> 
</p>

Overall, the OMP method seems to be the best. There are two slightly different OMP methods: [iss_OMP](https://github.com/jduffield65/iss/blob/a8c4d104274775dae67c2bace447476f40f0f355/@iss_OMP/iss_OMP.m) and [iss_OMP_ConstantBackground_WeightDotProduct](https://github.com/jduffield65/iss/blob/a8c4d104274775dae67c2bace447476f40f0f355/@iss_OMP_ConstantBackground_WeightDotProduct/iss_OMP_ConstantBackground_WeightDotProduct.m). The major differences are that for every iteration, iss_OMP, refits the coefficient for each background vector as well as for each gene already added, whereas iss_OMP_ConstantBackground_WeightDotProduct fits the background at the beginning and doesn't update it thereafter. Also, the score used to find the next best gene is different in iss_OMP_ConstantBackground_WeightDotProduct, and is explained in this document: [OMP Maths.md](https://github.com/jduffield65/iss/blob/a8c4d104274775dae67c2bace447476f40f0f355/@iss_OMP_ConstantBackground_WeightDotProduct/OMP%20Maths.md). These couple of tweaks seem to make iss_OMP_ConstantBackground_WeightDotProduct the better choice.

By default, bridge_process_template.m uses iss_OMP. To use iss_OMP_ConstantBackground_WeightDotProduct, change the [first line](https://github.com/jduffield65/iss/blob/4e0d03d53ad006c92073db3f20b9f6fd21557f0a/bridge_process_template.m#L8) from ```o = iss_OMP;``` to ```o = iss_OMP_ConstantBackground_WeightDotProduct;```. The default thresholds for this method are also different. 

### Viewing specific genes
To see the distribution of a specific gene or specific set of genes, run ```iss_change_plot(o,CallSpotsMethod,'Neuron',GeneNames)``` with the plot open, where ```CallSpotsMethod``` is ```'Prob'```, ```'DotProduct'```, ```'Pixel'``` or ```'OMP'``` as before. 'Neuron' is there just to specify that the genes are of neuron type. Some CodeFiles also contain non-neurons in which case, changing this argument to 'NonNeuron' will show different genes. GeneNames is a cell array containing the names of the genes of interest, so to see Plp1 and somatostatin  with the probability method, run ```iss_change_plot(o,'Prob','Neuron',[{'Plp1'},{'Sst'}])```. The result is shown below.

<p float="left">
<img src="DebugImages/README/SstPlp1.png" width = "450"> 
</p>

The gene names given must exactly match those names in [```o.GeneNames```](https://github.com/jduffield65/iss/blob/2ec0f5fb924c28f76b06d5d9d00bc14f88d4b2ba/%40iss/iss.m#L541) which come from the codebook. To revert to showing all genes, run with ```GeneNames=o.GeneNames``` i.e. ```iss_change_plot(o,'Prob','Neuron',o.GeneNames)```. To see all genes except for Plp1 and somatostatin, run with ```GeneNames=setdiff(o.GeneNames,[{'Plp1'},{'Sst'}])```.

### Viewing specific spots
To see the distribution of a specific set of spots run ```iss_change_plot(o,CallSpotsMethod, 'Neuron', GeneNames, SpotSet)``` with the plot open. ```SpotSet``` is logical array and only spots ```s``` for which ```SpotSet(s) = 1``` are shown. This allows you to choose your own thresholding methods, which may differ from the [default ones](https://github.com/jduffield65/iss/blob/PixelBased/@iss/quality_threshold.m). An example dataset for which ```SpotSet = o.pxSpotScore>30 & o.pxSpotIntensity > 500;``` is shown below:

<p float="left">
<img src="DebugImages/README/SpecificSpots1.png" width = "450"> 
</p>

#### Clustered spots
You can also restrict the display to spots that are clustered, this acts as a guide to where the cell locations are. To do this, run ```SpotSetClustered = get_gene_clusters(o,CallSpotsMethod,r,k,SpotSet)``` followed by ```iss_change_plot(o,CallSpotsMethod,'Neuron',GeneNames, SpotSetClustered)```. A cluster is required to have ```k``` spots from ```SpotSet``` to be within a distance ```r``` pixels of each other. An example with ```CallSpotsMethod = Pixel```, ```r = 7```, ```k = 3``` and ```SpotSet = o.pxSpotScore>30 & o.pxSpotIntensity > 500;``` is shown below:

<p float="left">
<img src="DebugImages/README/SpecificSpots2.png" width = "450"> 
</p>

By default (running ```get_gene_clusters(o,'Pixel')```), ```r = 18```, ```k = 20``` and ```SpotSet = o.quality_threshold(CallSpotsMethod)```. This is shown below for the pixel based method:

<p float="left">
<img src="DebugImages/README/SpecificSpots3.png" width = "450"> 
</p>

### Visualising individual spots
To view the dot product assignment of a particular gene, with the plot open, run [```iss_view_codes(o,234321,Norm)```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L119). This will cause a crosshair to appear on the plot, then just click on the spot of interest as shown below.

<p float="left">
<img src="DebugImages/README/Crosshair.png" width = "450"> 
</p>

234321 is just the figure number of the plot (should always be the same). Norm controls the normalisation applied to the spot and gene codes. You can set Norm equal to 1,2 or 3 to highlight certain features:
* Norm = 1: This gives the raw values. For the spot in the previous section, this would be:

<p float="left">
<img src="DebugImages/README/Norm1.png" width = "450"> 
</p>

* Norm = 2: This normalises in the same way used in ```call_spots``` i.e. it [normalises by the percentile](https://github.com/jduffield65/iss/blob/3f0ec254b1bd71b1c4b15ebcf9b319e3fd82f70d/%40iss/call_spots.m#L42) given by [```o.SpotNormPrctile```](https://github.com/jduffield65/iss/blob/59a7583fef8bd0231cbc0182394fcdcff0c84a9c/%40iss/iss.m#L267) in each colour channel and round and then depending on ```o.CallSpotsCodeNorm```, it [normalises the resultant so the either the whole code has L2 norm of 1 or each round does](https://github.com/jduffield65/iss/blob/3f0ec254b1bd71b1c4b15ebcf9b319e3fd82f70d/%40iss/call_spots.m#L163-L187). The plot in the previous section used Norm = 2.

* Norm = 3: This normalises [each colour channel by the percentile given by ```o.SpotNormPrctile```  across all rounds](https://github.com/jduffield65/iss/blob/3f0ec254b1bd71b1c4b15ebcf9b319e3fd82f70d/iss_view_prob.m#L41-L49). Using this, the spot we are considering would appear like this:

<p float="left">
<img src="DebugImages/README/Norm3.png" width = "450"> 
</p>

Also, if you know the index (```SpotNo```) of the spot you are interested in but don't want to find it in the plot, you can just run ```iss_view_codes(o,234321,Norm,SpotNo)```.

The equivalent function for the probability method is [```iss_view_prob(o,234321,Norm,CallSpotsMethod)```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L123) which is called in exactly the same way, only can specify whether ```CallSpotsMethod``` is ```'Prob'``` or ```'Pixel'```. If this is not specified, it will default to the method used in the plot currently open. The figure is the same except the gene shown is that given by ```o.pSpotCodeNo``` or  ```o.pxSpotCodeNo``` not ```o.SpotCodeNo```. Also, there is an extra plot which shows you how the overall probability is made up from the individual probabilities of each round and channel, relative to the probability that the spot can be explained by the background alone (without any genes). This plot for the example spot is given below:

<p float="left">
<img src="DebugImages/README/view_prob.png" width = "550"> 
</p>

The plot also gives the spot score, log probability over background, score deviation and intensity values for its assignment. Log probability over background is the total probability relative to the probability given by just the background distribution i.e. sum of all squares in bottom plot. Score deviation is the standard deviation of the log probabilities for assigning the spot to every gene in the codebook. This is included as if ```o.pSpotScoreDev+o.pSpotScore<o.pDevThresh```, then the assignment is rejected. This is to get rid of spots that have a similar probability when matched to every gene. Also, the values of these variables are coloured - green means that value caused the gene to be automatically accepted, red means that that value caused the gene to fail (green supersedes red most of the time). If all variables are black, then the match is also accepted.

### Understanding the probability method
If you left/right click on a particular square in the bottom plot of the ```iss_view_prob``` window, the plot on the left/right below will appear, further explaining the probability for that particular round and channel. These particular plots are for round 7, channel 2 for the above example spot. These are the absolute probabilities (i.e. without subtracting the probability due to background alone).

<p float="left">
<img src="DebugImages/README/NonCodeGoodPlt1.png" width = "400"> 
<img src="DebugImages/README/NonCodeGoodPlt2.png" width = "400">
</p>

The left plot shows how the probability for the round/channel that was clicked varies with spot intensity. The red line is the intensity for the spot considered, so where this intersects the blue curve indicates the absolute probability of the square that was clicked (-3.284 here). This then basically tells you how well the assign gene (blue curve) explains the spot intensity (red line) in this particular round/channel, the closer the red line is to the peak of the blue curve, the better.

The right plot reveals a bit more about what is happening behind the curtain. It shows the underlying probability distributions that lead to the final probability. The probability that we are trying to maximise is the probability that the spot, <img src="https://tex.s2cms.ru/svg/s" alt="s" />, can be explained by the gene, <img src="https://tex.s2cms.ru/svg/g" alt="g" />: <img src="https://tex.s2cms.ru/svg/P(s%5Cmid%20g)" alt="P(s\mid g)" />.

The equation 

<img src="https://tex.s2cms.ru/svg/P(s%5Cmid%20g)%20%3D%20%5Cint%20P(%5Clambda)P_b(s-%5Clambda%20g)d%5Clambda" alt="P(s\mid g) = \int P(\lambda)P_b(s-\lambda g)d\lambda" />

is saying that the probability of matching spot <img src="https://tex.s2cms.ru/svg/s" alt="s" /> to gene <img src="https://tex.s2cms.ru/svg/g" alt="g" /> is equal to the probability that the intensity of <img src="https://tex.s2cms.ru/svg/s" alt="s" /> once the intensity of the gene scaled by <img src="https://tex.s2cms.ru/svg/%5Clambda" alt="\lambda" /> has been removed can be explained by the background . This is then weighted by the probability of that scaling, <img src="https://tex.s2cms.ru/svg/P(%5Clambda)" alt="P(\lambda)" /> and then summed over all possible scalings, <img src="https://tex.s2cms.ru/svg/%5Clambda" alt="\lambda" />.

The probability distribution of <img src="https://tex.s2cms.ru/svg/%5Clambda" alt="\lambda" /> is different depending whether the particular round/channel, <img src="https://tex.s2cms.ru/svg/%5Br%2Cb%5D" alt="[r,b]" />, appears in ```o.CharCodes(g)``` or not, for each gene, there is one such colour channel for each round so ```o.nRounds``` incidences in total. We call this set <img src="https://tex.s2cms.ru/svg/S_g" alt="S_g" />, it is different for each gene <img src="https://tex.s2cms.ru/svg/g" alt="g" />. For <img src="https://tex.s2cms.ru/svg/%5Br%2Cb%5D" alt="[r,b]" /> not in <img src="https://tex.s2cms.ru/svg/S_g" alt="S_g" />, we expect the spot not to appear so the probability is peaked at <img src="https://tex.s2cms.ru/svg/%5Clambda%3D0" alt="\lambda=0" />. The actual distribution we use is (<img src="https://tex.s2cms.ru/svg/C_1" alt="C_1" /> is a constant):

<img src="https://tex.s2cms.ru/svg/P_%7B%5Br%2Cb%5D%5Cnot%5Cin%20S_g%7D(%5Clambda)%20%3D%20%5Cfrac%7BC_1%7D%7B2%7De%5E%7B-C_1%5Cmid%5Clambda%5Cmid%7D" alt="P_{[r,b]\not\in S_g}(\lambda) = \frac{C_1}{2}e^{-C_1\mid\lambda\mid}" />

The round 7, channel 2 example shown, corresponds to this set, so the blue curve on the right follows this distribution.

On the other hand, for <img src="https://tex.s2cms.ru/svg/%5Br%2Cb%5D" alt="[r,b]" /> that are in <img src="https://tex.s2cms.ru/svg/S_g" alt="S_g" />, we do expect the spot to appear and thus the scaling of the gene cannot be 0 or negative. We use the rayleigh distribution, which is such that <img src="https://tex.s2cms.ru/svg/P(%5Clambda%5Cleq%200)%20%3D%200" alt="P(\lambda\leq 0) = 0" /> and also the probability stays near the peak level for quite a large range of <img src="https://tex.s2cms.ru/svg/%5Clambda" alt="\lambda" /> which is desirable as to not penalise spot intensities that are higher than expected too much. The actual distribution we use is (<img src="https://tex.s2cms.ru/svg/C_2" alt="C_2" /> is a constant):

<img src="https://tex.s2cms.ru/svg/P_%7B%5Br%2Cb%5D%5Cin%20S_g%7D(%5Clambda%5Cleq%200)%20%3D%200%20%5C%5C%20P_%7B%5Br%2Cb%5D%5Cin%20S_g%7D(%5Clambda%20%3E0)%20%3D%20%5Cfrac%7B%5Clambda%7D%7BC_2%5E2%7De%5E%7B-%5Clambda%5E2%2F(2C_2%5E2)%7D" alt="P_{[r,b]\in S_g}(\lambda\leq 0) = 0 \\ P_{[r,b]\in S_g}(\lambda &gt;0) = \frac{\lambda}{C_2^2}e^{-\lambda^2/(2C_2^2)}" />

For the example spot, ```[round 2, colour channel 2]``` is in the set <img src="https://tex.s2cms.ru/svg/S_g" alt="S_g" /> and the plots are shown below:

<p float="left">
<img src="DebugImages/README/CodeGoodPlt1.png" width = "400"> 
<img src="DebugImages/README/CodeGoodPlt2.png" width = "400">
</p>

The constants <img src="https://tex.s2cms.ru/svg/C_1" alt="C_1" /> and <img src="https://tex.s2cms.ru/svg/C_2" alt="C_2" /> are [```o.ExpConst```](https://github.com/jduffield65/iss/blob/d8dc313f1a50d47a1df386098bcea2811f09dbf6/%40iss/iss.m#L607) and [```o.RaylConst```](https://github.com/jduffield65/iss/blob/d8dc313f1a50d47a1df386098bcea2811f09dbf6/%40iss/iss.m#L603) respectively. The default values were chosen so the mean values of the probability distributions matched the mean values of <img src="https://tex.s2cms.ru/svg/%5Clambda" alt="\lambda" /> in a particular data set.

The red curve on the right plot represents the background distribution which is just the [histogram of the raw data after the initial filtering](https://github.com/jduffield65/iss/blob/ffcffdaa492369e14ea4cfef214025e84e1becdf/%40iss/extract_and_filter.m#L152). This tends to be strongly peaked where the spot intensity is equal to the scaled intensity. Because the histogram only takes discrete integer values, we need to turn the integral into a sum, and also making the substitution <img src="https://tex.s2cms.ru/svg/x%3D%5Clambda%20g" alt="x=\lambda g" />, we get (The <img src="https://tex.s2cms.ru/svg/1%2Fg" alt="1/g" /> factor clearly blows up for <img src="https://tex.s2cms.ru/svg/g%3C1" alt="g&lt;1" /> so we also normalise <img src="https://tex.s2cms.ru/svg/%5Cfrac%7B1%7D%7Bg%7DP%5Cleft(%5Cfrac%7Bx%7D%7Bg%7D%5Cright)" alt="\frac{1}{g}P\left(\frac{x}{g}\right)" /> so it has a sum of 1 over all <img src="https://tex.s2cms.ru/svg/x" alt="x" />):

<img src="https://tex.s2cms.ru/svg/P(s%5Cmid%20g)%20%3D%20%5Cfrac%7B1%7D%7Bg%7D%5Csum_%7Bx%7D%20P%5Cleft(%5Cfrac%7Bx%7D%7Bg%7D%5Cright)P_b(s-x)" alt="P(s\mid g) = \frac{1}{g}\sum_{x} P\left(\frac{x}{g}\right)P_b(s-x)" />

Thus, from this right hand plot, to get the absolute probability of the clicked upon square, we need to multiply the two curves together and then sum the resultant over all x. This means the degree of overlap determines the probability, the greater the overlap the better. For example, looking at the ```[round 7, colour channel 2]``` case, the only way the log probability can get closer to the peak in the left hand plot is by the spot intensity reducing a little. The effect of this on the right hand plot would be to shift the red curve to the left which clearly increases the degree of overlap with the blue curve. Both the rounds/channels considered so far show pretty good overlap, an example with considerably worse overlap is ```[round 7, colour channel 3]``` shown below:

<p float="left">
<img src="DebugImages/README/NonCodeBadPlt1.png" width = "400"> 
<img src="DebugImages/README/NonCodeBadPlt2.png" width = "400">
</p>

This example also exhibts a potential pitfall of the method in dealing with bleedthrough. ```[round 2, colour channel 3]``` is not in
<img src="https://tex.s2cms.ru/svg/S_g" alt="S_g" /> for Nrn1 so we assume <img src="https://tex.s2cms.ru/svg/%5Clambda" alt="\lambda" /> is most likely to be 0, this means for a spot to match to Nrn1, we would expect that spot to have intensity of zero in ```[round 2, colour channel 3]```. But looking at the predicted code for Nrn1, the gene intensity in ```[round 2, colour channel 3]``` is 1891, so why would be possibly expect the spot to have zero intensity? This seems to throw away information, we have learned about the bleedthrough between colour channel 2 and 3. The Bleedthrough does have some effect though. Because of the <img src="https://tex.s2cms.ru/svg/1%2Fg" alt="1/g" /> factor, the <img src="https://tex.s2cms.ru/svg/%5Cfrac%7B1%7D%7Bg%7DP%5Cleft(%5Cfrac%7Bx%7D%7Bg%7D%5Cright)" alt="\frac{1}{g}P\left(\frac{x}{g}\right)" /> curve is flattened out. This reduces the peak probability and the magnitude of the gradient (to get from the peak log probability to the peak probability - 1 for ```[round 2, colour channel 3]``` would require a change in spot intensity of around 700 whereas in ```[round 7, colour channel 2]``` where there is no bleedthrough, this would require a change of 8). From this we see that the probability of rounds/channels with high bleedthrough are insensitive to spot intensity but have a low probability (although with bleedthrough, the background distribution tends to have a wider peak which helps nullify the reduction in peak probability - the wider the background distribution peak, the more overlap with the wide <img src="https://tex.s2cms.ru/svg/%5Cfrac%7B1%7D%7Bg%7DP%5Cleft(%5Cfrac%7Bx%7D%7Bg%7D%5Cright)" alt="\frac{1}{g}P\left(\frac{x}{g}\right)" /> caused by bleedthrough, thus the higher the max probability). To counteract this, we could add these rounds/channels with high bleedthrough to <img src="https://tex.s2cms.ru/svg/S_g" alt="S_g" />, but this would make the size of <img src="https://tex.s2cms.ru/svg/S_g" alt="S_g" /> vary between genes which may introduce a subtelty when comparing the probabilities of assigning to genes, i.e. an artificial spot with high postive intensity in all rounds/channels would preferentially match with genes for which the size of <img src="https://tex.s2cms.ru/svg/S_g" alt="S_g" /> is larger.

### Visualising filtering step
The tiles are filtered in the [```extract_and_filter``` step](https://github.com/jduffield65/iss/blob/34786306a0ee42fcea10188b7c5c688362eb20b6/%40iss/extract_and_filter.m#L230-L239). The goal of this step is to emphasize the spots over the background. To see if it has worked as intended or if the filtering parameters are correct, you can run [```view_filtering(o,r,t)```](https://github.com/jduffield65/iss/blob/PixelBased/view_filtering.m); r is the round of interest and t is the tile of interest. You can run this before the tiles have been produced, the only requirement is that the ```o``` object must have the following properties specified:
* ```o.InputDirectory```
* ```o.FileBase```
* ```o.RawFileExtension```
* ```o.AnchorRound```
* ```o.AnchorChannel```
* ```o.DapiChannel```
* ```o.TileSz```

The first image that will appear is the raw image for the first colour channel. You can use the horizontal scroll bar to change colour channel, e.g. an example anchor image is shown below on the left. You can then press the Filter button to see the filtered image (i.e. what the files in ```o.TileDirectory``` will end up looking like).

<p float="left">
<img src="DebugImages/README/UnFiltered.png" width = "400"> 
<img src="DebugImages/README/Filtered.png" width = "400"> 
</p>

With the filter button pressed, you can change the radius of the filter used with the vertical slider and the plot should update automatically. The default value is the value of ```o.ExtractR1``` that would be used in the pipeline automatically. If you decide that another value is more suitable, you can just run in the command window: ```o.ExtractR1 = NewValue;```. The plots below show the same region that is unfiltered, filtered with filter radius of 3 (default) and 8 respectively.

<p float="left">
<img src="DebugImages/README/UnFilteredClose.png" width = "400"> 
<img src="DebugImages/README/FilteredCloseR3.png" width = "400"> 
<img src="DebugImages/README/FilteredCloseR8.png" width = "400"> 
</p>


If you are viewing the Dapi colour channel in the Anchor round, then the vertical slider controls ```o.DapiR``` instead and the [filtering is different](https://github.com/jduffield65/iss/blob/34786306a0ee42fcea10188b7c5c688362eb20b6/%40iss/extract_and_filter.m#L233-L234). 

### OMP Visualisation
There are two functions to help debug omp gene assignments: [iss_view_omp](https://github.com/jduffield65/iss/blob/1e2a7e610b14e98d1859000396718c456c6324ef/@iss_OMP/iss_view_omp.m) and [iss_view_spot_omp2](https://github.com/jduffield65/iss/blob/1e2a7e610b14e98d1859000396718c456c6324ef/@iss_OMP/iss_view_spot_omp2.m). 

The first is run through ```iss_view_omp(o, 234321, Norm, IncludeGT, SpotNo)``` with the ```o.plot()``` window open. ```Norm``` is the same as with ```iss_view_codes```, ```IncludeGT``` is set to false by default and should remain as this. If you know the index of the spot within ```o.ompSpotGlobalYX``` that you are interested in or indeed the YX location of interest, you can set ```SpotNo``` to this value. Otherwise, leave this argument empty and the crosshair will appear on the plot.

An example with ```Norm = 2``` is shown below:

<p float="left">
<img src="DebugImages/README/ompView.png" width = "650"> 
</p>

The top plot just shows the spot color after the specified normalisation. The bottom plot shows the coefficients found by the OMP method for each of the genes as well as the background for this pixel. The second plot shows the predicted code which is the sum of all these coefficients multiplied by the bled codes. The third plot shows the breakdown of the score i.e. ompSpotScore is the sum of all squares in this image. For more details about ompScore look [here](https://github.com/jduffield65/iss/blob/1e2a7e610b14e98d1859000396718c456c6324ef/@iss_OMP_ConstantBackground_WeightDotProduct/Thresholding.md). If you left click on a gene in the bottom plot, the third plot will show the score due to that gene. If you right click, the second plot will change to the predicted code without that gene and the third plot will show the error (difference between top two plots) without that gene.

The OMP algorithm finds the coefficient of every gene at every pixel of the image, from which you can produce an image for each gene. ```iss_view_spot_omp2``` shows these coefficient images. It is run through ```iss_view_spot_omp2(o, 234321, ImSz, SpotLocation, ScoreMethod, Track, SpotNo)``` with the ```o.plot()``` window open. ```ImSz``` is the radius of the image produced. ```SpotLocation``` is true if you want the image to be centered on the closest spot to the crosshair and false if you want the exact position of the crosshair. ```ScoreMethod``` is the same as ```CallSpotsMethod``` in ```iss_change_plot```. It is only used to find the location, most likely it will be ```'OMP'```. ```Track``` can be set to true to give some extra plots if the class is ```iss_OMP``` but cannot be used if the class is ```iss_OMP_ConstantBackground_WeightDotProduct```. ```SpotNo``` is the same as with ```iss_view_omp```.

The plot for the same spot shown above is shown here:

<p float="left">
<img src="DebugImages/README/ompCoefImage.png" width = "650"> 
</p>

Red is a positive coefficient, blue is negative. A green circle is a local maxima, a grey ```x``` is a spot found by the omp algorithm and saved to ```o.ompSpotGlobalYX``` but failed the thresholding. A grey ```+``` is a saved spot that passed the thresholding. 


## View round/channel images of spot
Another function that is useful for evaluating the gene calling results is [```iss_view_spot```](https://github.com/jduffield65/iss/blob/e191656ae39590ceb3a7c9ea9da67ae3fb71b6fe/@iss_Base/iss_view_spot.m), which plots the location around a spot in each round and colour channel. It can be used in collaboration with any of the gene calling methods. It is run through: ```iss_view_spot(o, 234321, ImSz, SpotLocation, ScoreMethod, IncludeGT, Filter, Norm, SpotNo)```. Here, ```Norm``` is either ```true``` or ```false``` to plot the channel/round normalised or raw colors. ```Filter``` is either ```true``` or ```false``` to plot the colors after or before the ```extract_and_filter``` step. All the other arguments are as described in previous functions.

It is particularly useful to use and compare with ```iss_view_spot_omp2```, and the corresponding plot is shown below with ```Norm = false``` and ```Filter = true```:

<p float="left">
<img src="DebugImages/README/ompViewSpot.png" width = "950"> 
</p>

By comparing the two images, you can see the channel 0 background spot in the south-east quadrant as well as the background spot towards the centre in channels 1, 2, 4 and 5. The spot shape in these channels appears to be the same in all rounds, indicating they are background not genes. 

The Aldoc gene can be seen along the diagonal, in the plots with the green crosshair. It is quite hard to evaluate due to all the background but there does seem to be a greater intensity in the middle of these green channels in each round compared to the centre of other channels of the same round. In particular, round 7, channel 6 is quite convincing evidence for Aldoc as there is no background for channel 6 so the spot in this channel must be explained by a gene. 

You can see the Cadps2 (```CharCode = 2345601```) spot emerging in the south-east especially well in (round 2, channel 3); (round 5, channel 6) and (round 7, channel 1). The only other signal of note in the ```iss_view_spot_omp2``` plot is Gda (```CharCode = 1526304```) in the south-east quadrant. The signal in (round 4, channel 6) and maybe (round 5, channel 3) seem to differ from that expected from background and have a spot in the location of the supposed Gda. In all the other rounds though, either there is no spot or it coincides with background. Thus it is hard to determine if this spot is real. As can be seen from the grey ```x``` in the ```iss_view_spot_omp2``` plot, this spot did not pass the thresolding process, which I think is fair enough.

### ```iss_view_spot_omp3```
The functions ```iss_view_spot_omp2``` and ```iss_view_spot``` are explicitly combined in the function [```iss_view_spot_omp3```](https://github.com/jduffield65/iss/blob/c00c764053f6230adb44eb242ffc8bb4317aadff/@iss_OMP/iss_view_spot_omp3.m).  It is run through: ```iss_view_spot_omp3(o, 234321, ImSz, SpotLocation, ScoreMethod, Track, SpotNo, Norm)``` and both the two above plots will open. The arguments are all the same as for ```iss_view_spot_omp2``` except for the additional argument ```Norm``` which can be ```true``` or ```false``` as with ```iss_view_spot```. The default value is ```true``` for this function as it shows the effect of removing genes more clearly. The SpotColor image with ```Norm=true``` is shown below:

<p float="left">
<img src="DebugImages/README/ompViewSpot3_Original.png" width = "950"> 
</p>

If you then left click on a gene image in the coefficient plot that also opened, then the SpotColor image will adjust to remove the gene you clicked on. Clicking on any background gene will remove all background genes. Right clicking will reinstate the gene. The SpotColor without background and then without background and Aldoc are shown below:

<p float="left">
<img src="DebugImages/README/ompViewSpot3_RemoveBackground.png" width = "950"> 
</p>

<p float="left">
<img src="DebugImages/README/ompViewSpot3_RemoveBackgroundAldoc.png" width = "950"> 
</p>

So in a perfect example, each of the images will go completely white (to 0) after we have clicked on all the genes that are actually there. 
