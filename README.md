# Robust_ranking
This is the example code of [1]Fu et al. Robust Subjective Visual
Property Prediction from Crowdsourced Pairwise Labels, IEEE TPAMI, to
appear.
[2] Fu et al. Interestingness prediction by robust learning to rank,
ECCV 2014.

I have the following codes
1) "odet.R" use Cross-validation to select lambda in solving HLASSO.
2) "fyw.R" simulates 16 nodes with 3000 random comparisons with correct
direction then reverse 10% of them. Binary choice is used for each
comparison, which may be similar to your case. HLASSO is used to detect
outliers. Here CV doesn't pick out any outlier. But you can still
delete top p% of outliers.
There are more details in the annotation of codes.
3) our codes on scene and pubfig. (I dont have these codes in this
machine; but luckly, I searched my gmail and found I previous sent the
similar codes to Jiechao for help. I use that version.)

4) Age_exp_v3.R is the version of using low-level feature for pruning.
Thanks!


Please let me know if there is any other problem. 

Yanwei Fu (y.fu(at)qmul dot ac dot uk;) 
ztwztq2006(at) gmail dot com
