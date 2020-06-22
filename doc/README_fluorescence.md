# easyCLIP: adapter fluorescence calculations

From the notebook nonclip/adapter_quant/STDs experiments.ipynb:

To show a linear function of signal vs fmols in staples:

```python
import sameRiver
import sameRiver.adapter_fluorescence as af

params5, params3 = af.get_parameters_and_draw_staple(fname='180523.xlsx', sheet_name='180523_1')
```

The outputs of this are staple_aL5_signal_vs_fmols.pdf and staple_aL3_signal_vs_fmols.pdf.
These are fits of fluorescence in the staple oligos to a line.


From the notebook nonclip/adapter_quant/STDs experiments.ipynb:

For the accuracy of the fit and parameters on the sample used to get parameters:

```python
import sameRiver
import sameRiver.adapter_fluorescence as af

# File with the samples used for parameters.
df = af.load_sheet(fname='180522.xlsx', sheet_name='180522_1')

df = df[df['fmols']<100]
df = df[df['fmols']>3]

# Get fit and parameters
df = af.est_fmols_linear(df, staples=None)

l5 = df[df['Object']=='αL5']
l3 = df[df['Object']=='αL3']

# This fit is displayed with these calls:
plt.plot(
    l5['fmols'].tolist(), l5['Est. fmols'].tolist(), 'r.', markersize=14, alpha=0.3,
    markeredgewidth=0)
plt.plot(
    l3['fmols'].tolist(), l3['Est. fmols'].tolist(), 'g.', markersize=14, alpha=0.3,
    markeredgewidth=0)

# Errors are obtained from:

def errors(df):
    df['% error'] = [abs(100*(est-fmols)/fmols) for est, fmols in zip(
         df['Est. fmols'], df['fmols'])]
    df['fmols error'] = [abs(est-fmols) for est, fmols in zip(
         df['Est. fmols'], df['fmols'])]
    return df

l5 = errors(l5)
l3 = errors(l3)
```

For the usage of the adapter shift method for quantification, code is in
 nonclip/adapter_quant/L5,L3 by aL5,aL3 estimattion.ipynb.

Specifically, we get the Ligation efficiencies:
```python
import sameRiver
import sameRiver.adapter_fluorescence as af
import importlib
importlib.reload(af)

obj_color = af.obj_color

exp_folder = '/Users/dfporter/pma/dataAndScripts/clip/experiments/exp35 hnRNPC FBL AURKA RPS3/'

# Get replicate 3 from 180525. Not used for some reason.
df_r3 = af.load_sheet(
    fname=exp_folder+'/exp35.xlsx',
#    fname='../Shift based quantification protocol and results.xlsx',
    sheet_name='180525 hnRNPC R3 qRNA')

# Get replicate 3. Not used for some reason.
df2 = af.load_sheet(
    fname=exp_folder+'/exp35.xlsx',
#    fname='../Shift based quantification protocol and results.xlsx',
    sheet_name='qRNA hnRNPC R3'
    #sheet_name='qRNA hnRNPC R1,R2'
)

# Get replicates 1,2,3 from 180501. This is all the data.
df3 = af.load_sheet(
    fname=exp_folder+'/exp35.xlsx',
    sheet_name='180501 hnRNPC R1,2,3 qRNA')
#df = pandas.concat([df2, df3])

print('&' * 14)

def ests_Kas(df):
    [[m5, b5], [m3, b3]] = af.get_parameters(fname='180522.xlsx', sheet_name='180522_1')
    df['Replicate'] = ['R' + str(x)[:1] for x in df.Replicate]
    # The staples parameter is very imporant. The ligation efficiencies are +/- ~20%
    # based on whether staples=None. If staples=None, then fluorescence values are
    # compared to the run used for parameters (180522_1) without normalization.
    # If staples is not None, then values are scaled based on relative fluorescence
    # of the antisense oligos in the 50 fmol staples lanes. It is not clear which
    # is the better method, but we have opted to use scaling, as it makes perhaps a little more sense.
    af.est_fmols_linear_using_params(df, [[m5, b5], [m3, b3]], staples=True)
    
    return df

df3 = ests_Kas(df3)

eff = ligation_efficiency(df3)
print("Lig efficiencies:")
print(eff)
```

Ligation efficiencies are then written in by hand to nonclip/adapter_quant/all_ligation_eff_estimates.xlsx.

The final ligation efficiency graph is output by nonclip/adapter_quant/Ligation efficiency.ipynb,
 based on nonclip/adapter_quant/all_ligation_eff_estimates.xlsx.

