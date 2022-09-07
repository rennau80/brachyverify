#!/usr/bin/python
import pylab as plt
import pydicom as dcm
import numpy as np
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import sys 
import os
import datetime

matplotlib.rc('xtick', labelsize=7) 
matplotlib.rc('ytick', labelsize=7) 

def modification_date(filename):
    t = os.path.getmtime(filename)
    return datetime.datetime.fromtimestamp(t)

def calc_dose_single_position(sourcez,sourcey,sourcex,stime,source_rad,z,y,x):
    dist=np.sqrt((sourcez-z)**2+(sourcey-y)**2+(sourcex-x)**2)
    if dist<=1:
        return stime*source_rad
    else:
        return stime*source_rad/(dist*dist)

def calc_source_angle_and_radius(xbefore,ybefore,zbefore,xnext,ynext,znext,sourcex,sourcey,sourcez,x,y,z):
    source_dx=xbefore-xnext
    source_dy=ybefore-ynext
    source_dz=zbefore-znext
    radial_dx=x-sourcex
    radial_dy=y-sourcey
    radial_dz=z-sourcez

    radius=np.sqrt(radial_dx**2+radial_dy**2+radial_dz**2)

    vect_abs=source_dx*radial_dx+source_dy*radial_dy+source_dz*radial_dz
    vect1_abs=source_dx**2+source_dy**2+source_dz**2
    vect2_abs=radial_dx**2+radial_dy**2+radial_dz**2
    angle=np.degrees(np.arccos(vect_abs/(np.sqrt(vect1_abs*vect2_abs))))

    if radius>150:radius=150 #as dose only 15cm around in anisotrophy file, so take 15cm value for all r>150
    return radius,angle

def calc_anisotropie(radius,angle):
    f_anisotr=open('anisotropie.txt','r')
    line_dummy=0
    for l in f_anisotr:
        if line_dummy>0:
            l_=l.split('\t')
            if np.abs(radius-np.float(l_[0]))<2.5: angle_values=l_[1:-1]
        line_dummy+=1
    angle_index=int(angle/180.*len(angle_values))
    f_anisotr.close()
    return np.float(angle_values[angle_index])

def calc_dose_single_position_radial(xbefore,ybefore,zbefore,xnext,ynext,znext,sourcez,sourcey,sourcex,stime,source_rad,z,y,x):
    source_dx=xbefore-xnext
    source_dy=ybefore-ynext
    source_dz=zbefore-znext

    radial_dx=x-sourcex
    radial_dy=y-sourcey
    radial_dz=z-sourcez

    vect_abs=np.sqrt(source_dx*radial_dx+source_dy*radial_dy+source_dz*radial_dz)
    vect1_abs=np.sqrt(source_dx**2+source_dy**2+source_dz**2)
    vect2_abs=np.sqrt(radial_dx**2+radial_dy**2+radial_dz**2)
    
    angle=np.degrees(np.arctan(vect_abs/(vect1_abs*vect2_abs)))

    dist=np.sqrt((sourcez-z)**2+(sourcey-y)**2+(sourcex-x)**2)
    if dist<=2:
        return stime*source_rad
    else:
        return stime*source_rad/(dist*dist)

def extr_dose_info(datastring):
    pos1=str.find(datastring,'DS')
    pos2=str.find(datastring[pos1:],'\n')+pos1
    pos3=datastring[pos1:pos2].find('"')
    pos4=datastring[pos1+pos3+1:pos2].find('"')
    id_tmp=datastring[pos1+pos3+1:pos1+pos3+1+pos4]
    return id_tmp

def extr_contour_info(datastring):
    pos1=str.find(datastring,'DS')
    pos2=str.find(datastring[pos1:],'\n')+pos1
    pos3=datastring[pos1:pos2].find('[')
    pos4=datastring[pos1+pos3+1:pos2].find(']')
    id_tmp=datastring[pos1+pos3+1:pos1+pos3+1+pos4]
    return id_tmp

DEBUG=0

rp_file=sys.argv[-2]
rd_file=sys.argv[-1]

data_RP=dcm.read_file(rp_file)
data_rtdose=dcm.read_file(rd_file)

### get patient details ###
patname=str(data_RP.PatientName)
patname=patname.split('^')
patname=patname[0]+', '+patname[1]

patid=str(data_RP.PatientID)

#da es diesen eintrag nur bei bestaetigten plaenen gibt, nun nicht mehr drin
#patreviewdate=str(data_RP.ReviewDate)

data_RP_ASCII=str(data_RP)
start_kath_data=str.find(data_RP_ASCII,'300a, 0230')
fillval=-9999

#increase grid to side of min/max values for x,y,z, e.g. 10 meaning 5mm to each side
grid_increase=20
#resolution of dose grid in this calculation [mm]
grid_res=3

############################################
################# OTP DATA #################
############################################

print('read in OTP dose data ...')
print(patname,'[',patid,']')
x_dose=np.zeros(len(data_rtdose.RTDoseROISequence))
y_dose=np.zeros(len(data_rtdose.RTDoseROISequence))
z_dose=np.zeros(len(data_rtdose.RTDoseROISequence))
otp_dose_values=np.zeros(len(data_rtdose.RTDoseROISequence))
for curr_dose in range(len(data_rtdose.RTDoseROISequence)):
    tmp=extr_dose_info(str(data_rtdose.RTDoseROISequence[curr_dose]))
    otp_dose_values[curr_dose]=np.float(tmp)
    coordinates=extr_contour_info(str(data_rtdose.ROIContourSequence[curr_dose]))
    tmp=coordinates.split(',')
    x_dose[curr_dose]=np.float(tmp[0])
    y_dose[curr_dose]=np.float(tmp[1])
    z_dose[curr_dose]=np.float(tmp[2])

print('finished...')

#####################################################
########### SOURCE DATA #############################
#####################################################
#TO DO
#brachy appl setup dose and fractions planned
source_data1=data_RP.FractionGroupSequence[0]
#target prescr dose
source_data2=data_RP.DoseReferenceSequence
source_data3=data_RP.InstanceCreationDate
source_data4=data_RP.InstanceCreationTime

#sourceData
source_data5=data_RP.SourceSequence[0]

source_data1=str(source_data1)
source_data2=str(source_data2)
source_data3=str(source_data3)
source_data4=str(source_data4)
source_data5=str(source_data5)

source_data5=source_data5.split('\n')
source_output_data='\n'
for dummy in range(len(source_data5)):
    tmp=source_data5[dummy]
    tmp=tmp.split(' ')
    source_output_data=source_output_data+tmp[2]+' '+tmp[3]+' '+tmp[-1]+'\n'
source_data5=source_output_data

### willkuerlicher wert, wird spaeter angepasst anhand berechneter dosisdifferenzen:
source_rad=1.

#get channel number
ch_nr=str.find(data_RP_ASCII,'300a, 0280')
ch_nr_str=data_RP_ASCII[ch_nr:ch_nr+100]
stage=0
CN=""
for dummy in range(len(ch_nr_str)):
    if ch_nr_str[dummy]=="C": stage=1
    if stage==1:
        if np.char.isnumeric(ch_nr_str[dummy]):
            CN=ch_nr_str[dummy:dummy+3]
            if np.char.isnumeric(CN[-1])==0:
                CN=ch_nr_str[dummy:dummy+2]
            CN=np.int(CN.strip())
            stage=2

print('\n# of channels:',CN)

#obtain catheter data:
start_idx_kath=[]
shifter=0
for cath in range(CN):
    start_idx_kath_tmp=np.int(str.find(data_RP_ASCII,'300a, 0282) Channel Number                      IS: \"'+str(cath+1)+'\"'))
    start_idx_kath.append(start_idx_kath_tmp-70) #-70 to include number of control Points
start_idx_kath.append(start_idx_kath[-1]+30000)

#obtain 3D position
cath_x=np.ones((CN,100))*fillval
cath_y=np.ones((CN,100))*fillval
cath_z=np.ones((CN,100))*fillval
cath_time=np.ones((CN,100))*fillval
cath_total_time=np.ones((CN))*fillval

max_rad_time=0
max_rad_time_cath=''
rad_time_positions=[]
for cath in range(CN):
    cath_data=data_RP_ASCII[start_idx_kath[cath]:start_idx_kath[cath+1]]
    stage=1
    contr_point=0
    shifter=0
    #final cum time weight
    tmp1=cath_data.find('Channel Total Time')
    tmp2=cath_data[tmp1:tmp1+200]
    tmp2=tmp2[:tmp2.find('\n')-1]

    #####TEST wenn ein katheter komplett doppelte Zeit hat:
    err_time=1.

    cath_total_time[cath]=err_time*np.float(tmp2[tmp2.find('\"')+1:])
    while stage==1:
        if str.find(cath_data[shifter:],'(300a, 0112) Control Point Index                 IS: \"'+str(contr_point)+'\"')>0:
            tmp1=str.find(cath_data,'(300a, 0112) Control Point Index                 IS: \"'+str(contr_point)+'\"') 
            data3d=cath_data[tmp1:tmp1+310]
            # positions:
            tmp2=str.find(data3d,'3D')
            pos_data=data3d[tmp2:tmp2+70].strip()
            pos_data=pos_data[pos_data.find('['):]
            cath_x[cath,contr_point]=pos_data[1:pos_data.find(',')]
            cath_y[cath,contr_point]=pos_data[pos_data.find(',')+1:pos_data.find(',')+1+7]
            pos_data=pos_data.split(' ')
            cath_z[cath,contr_point]=pos_data[2][:7]

            # time weight
            tmp2=str.find(data3d,'Time Weight')
            time_data=data3d[tmp2:tmp2+50].strip()

            #TEST WIE STARK ZU SEHEN WENN EINZELNE POS FALSCH BERECHNET
            err_time=1
            #if cath==2 and contr_point==2: err_time=2

            #print(np.round(np.float(time_data[time_data.find('DS')+5:time_data.find('\n')-1]),2))
            if np.float(time_data[time_data.find('DS')+5:time_data.find('\n')-1])>max_rad_time:
                max_rad_time=np.float(time_data[time_data.find('DS')+5:time_data.find('\n')-1])
                max_rad_time_cath='Nr: '+str(cath+1)+', Pos: '+str(contr_point+1)
                cath_max_time_nr=cath
                cath_max_time_pos=contr_point
            cath_time[cath,contr_point]=err_time*np.float(time_data[time_data.find('DS')+5:time_data.find('\n')-1])

            if cath_time[cath,contr_point]>0.01: rad_time_positions.append(np.float(time_data[time_data.find('DS')+5:time_data.find('\n')-1]))

            shifter=tmp1+10
            contr_point+=1
        else:
            stage=0

text_cath_time='\n'
rad_time=0
for dummy in range(CN):
    print('cath #'+str(dummy+1)+' total time: '+ str(cath_total_time[dummy]))
    rad_time+=np.float(cath_total_time[dummy])
    text_cath_time=text_cath_time+' #'+str(dummy+1)+': '+str(np.round(cath_total_time[dummy],1))+'s, '
    if dummy>0 and dummy%5==0: text_cath_time=text_cath_time+'       '+'\n'
  
#print min/max values of x,y,z
min_x=-1*fillval
min_y=-1*fillval
min_z=-1*fillval
max_x=fillval
max_y=fillval
max_z=fillval

cathposnr=[]
cathposnrrunner=0
cathcumtime=[]
for cath in range(CN):
    cathposnrrunner=0
    for dummy in range(100):
        if cath_x[cath,dummy]!=fillval: 
            cathposnrrunner+=1
            tmp=cath_time[cath,dummy]
            if cath_x[cath,dummy]<min_x:min_x=cath_x[cath,dummy]
            if cath_x[cath,dummy]>max_x:max_x=cath_x[cath,dummy]
            if cath_y[cath,dummy]<min_y:min_y=cath_y[cath,dummy]
            if cath_y[cath,dummy]>max_y:max_y=cath_y[cath,dummy]
            if cath_z[cath,dummy]<min_z:min_z=cath_z[cath,dummy]
            if cath_z[cath,dummy]>max_z:max_z=cath_z[cath,dummy]
    cathcumtime.append(tmp)
    cathposnr.append(cathposnrrunner)

print('cath positions: ',cathposnr)
print('MINMAX found in Plan Data (x),(y),(z): ('+str(min_x)+','+str(max_x)+'),('+
        str(min_y)+','+str(max_y)+'),('+str(min_z)+','+str(max_z)+')')
print('MINMAX found in Dose Data (x),(y),(z): ('+str(np.min(x_dose))+','+
        str(np.max(x_dose))+'),('+str(np.min(y_dose))+','+str(np.max(y_dose))+'),('+
        str(np.min(z_dose))+','+str(np.max(z_dose))+')')

###################################################
##############CALC DOSE 3D#########################
###################################################
# check for start and end values of x,y,z coordinates and create new empty dose
# matrix to fill with catheter point doses
x_range=np.round((max_x-min_x)+grid_increase,0)
y_range=np.round((max_y-min_y)+grid_increase,0)
z_range=np.round((max_z-min_z)+grid_increase,0)
xvals_new=np.arange(np.int(min_x)-np.int(grid_increase/2),min_x+x_range+np.int(grid_increase/2),grid_res)
yvals_new=np.arange(np.int(min_y)-np.int(grid_increase/2),min_y+y_range+np.int(grid_increase/2),grid_res)
zvals_new=np.arange(np.int(min_z)-np.int(grid_increase/2),min_z+z_range+np.int(grid_increase/2),grid_res)
nx=len(xvals_new)
ny=len(yvals_new)
nz=len(zvals_new)
print('array for dose calculation: ',nx,ny,nz)
dose_matrix_3D=np.zeros((nz,ny,nx))
dose_matrix_3d_all_points=[]
for cathnr in range(CN):
    print('...processing 3D-calculation of cath #'+str(cathnr+1))
    for cathpos in range(np.int(cathposnr[cathnr]/2)):
        sx=cath_x[cathnr,cathpos*2]
        sy=cath_y[cathnr,cathpos*2]
        sz=cath_z[cathnr,cathpos*2]
        if cathpos==0:
            sbefx=cath_x[cathnr,cathpos*2]
            sbefy=cath_y[cathnr,cathpos*2]
            sbefz=cath_z[cathnr,cathpos*2]
            snex=cath_x[cathnr,(cathpos+1)*2]
            sney=cath_y[cathnr,(cathpos+1)*2]
            snez=cath_z[cathnr,(cathpos+1)*2]
        elif cathpos==np.int(cathposnr[cathnr]/2):
            sbefx=cath_x[cathnr,(cathpos-1)*2]
            sbefy=cath_y[cathnr,(cathpos-1)*2]
            sbefz=cath_z[cathnr,(cathpos-1)*2]
            snex=cath_x[cathnr,(cathpos)*2]
            sney=cath_y[cathnr,(cathpos)*2]
            snez=cath_z[cathnr,(cathpos)*2]
        else:
            sbefx=cath_x[cathnr,(cathpos-1)*2]
            sbefy=cath_y[cathnr,(cathpos-1)*2]
            sbefz=cath_z[cathnr,(cathpos-1)*2]
            snex=cath_x[cathnr,(cathpos+1)*2]
            sney=cath_y[cathnr,(cathpos+1)*2]
            snez=cath_z[cathnr,(cathpos+1)*2]

        st=(cath_time[cathnr,cathpos*2+1]-cath_time[cathnr,cathpos*2])/cathcumtime[cathnr]*cath_total_time[cathnr]

        if DEBUG:
            cmdose=calc_dose_single_position(sz,sy,sx,st,source_rad,sz+10,sy,sx)
            print('   #'+str(cathnr)+' pos'+str(cathpos)+' time: '+str(st)[:4]+'s 1cm dose: '+str(cmdose)[:4]+'Gy')

        if st!=0.0:
            for z in range(nz):
                for y in range(ny):
                    for x in range(nx):
                        #correct dose by anisotropie
                        (rad,ang)=calc_source_angle_and_radius(sbefx,sbefy,sbefz,snex,sney,snez,sx,sy,sz,
                                                               xvals_new[x],yvals_new[y],zvals_new[z])
                        red_factor=calc_anisotropie(rad,ang)
                        dose_matrix_3D[z,y,x]=dose_matrix_3D[z,y,x]+red_factor*calc_dose_single_position(sz,sy,sx,
                                                            st,source_rad,zvals_new[z],yvals_new[y],xvals_new[x])
                        dose_matrix_3d_all_points.append(dose_matrix_3D[z,y,x])

##############################################
##################### CALC DOSE ##############
##############################################
dose_matrix=np.zeros(len(otp_dose_values))
#run loop with filling dose matrix using every source position (see other code)
stage=1
check_points=np.arange(0,len(otp_dose_values))

for cathnr in range(CN):
    print('...processing point dose calculations of cath #'+str(cathnr+1))
    for cathpos in range(np.int(cathposnr[cathnr]/2)):
        sx=cath_x[cathnr,cathpos*2]
        sy=cath_y[cathnr,cathpos*2]
        sz=cath_z[cathnr,cathpos*2]
        if cathpos==0: 
            sbefx=cath_x[cathnr,cathpos*2]
            sbefy=cath_y[cathnr,cathpos*2]
            sbefz=cath_z[cathnr,cathpos*2]
            snex=cath_x[cathnr,(cathpos+1)*2]
            sney=cath_y[cathnr,(cathpos+1)*2]
            snez=cath_z[cathnr,(cathpos+1)*2]
        elif cathpos==np.int(cathposnr[cathnr]/2):
            sbefx=cath_x[cathnr,(cathpos-1)*2]
            sbefy=cath_y[cathnr,(cathpos-1)*2]
            sbefz=cath_z[cathnr,(cathpos-1)*2]
            snex=cath_x[cathnr,(cathpos)*2]
            sney=cath_y[cathnr,(cathpos)*2]
            snez=cath_z[cathnr,(cathpos)*2]
        else:
            sbefx=cath_x[cathnr,(cathpos-1)*2]
            sbefy=cath_y[cathnr,(cathpos-1)*2]
            sbefz=cath_z[cathnr,(cathpos-1)*2]
            snex=cath_x[cathnr,(cathpos+1)*2]
            sney=cath_y[cathnr,(cathpos+1)*2]
            snez=cath_z[cathnr,(cathpos+1)*2]

        st=(cath_time[cathnr,cathpos*2+1]-cath_time[cathnr,cathpos*2])/cathcumtime[cathnr]*cath_total_time[cathnr]

        if DEBUG: 
            cmdose=calc_dose_single_position(sz,sy,sx,st,source_rad,sz+10,sy,sx)
            print('   #'+str(cathnr)+' pos'+str(cathpos)+' time: '+str(st)[:4]+'s 1cm dose: '+str(cmdose)[:4]+'Gy')

        if st!=0.0:
            if cath_max_time_nr==cathnr and cath_max_time_pos-1==cathpos*2:
                red_factor=calc_anisotropie(8,90)
                cath_max_dose_5mm=red_factor*calc_dose_single_position(sz,sy,sx,st,source_rad,sz,sy+8,sx)
            for dummy in range(len(check_points)): 
                #correct dose by anisotropie
                (rad,ang)=calc_source_angle_and_radius(sbefx,sbefy,sbefz,snex,sney,snez,sx,sy,sz,x_dose[dummy],y_dose[dummy],z_dose[dummy])
                red_factor=calc_anisotropie(rad,ang)
                dose_matrix[dummy]=dose_matrix[dummy]+red_factor*calc_dose_single_position(sz,sy,sx,st,source_rad,z_dose[dummy],y_dose[dummy],x_dose[dummy])

#####################
#### dose scaling ###
#####################
corr_factors=[]
for dummy in range(len(check_points)):
    if otp_dose_values[dummy]<6:corr_factors.append(otp_dose_values[dummy]/dose_matrix[dummy])
corr_factor=np.average(corr_factors)    

for dummy in range(len(dose_matrix)):
    dose_matrix[dummy]=corr_factor*dose_matrix[dummy]

###############################
####extr Dosispunkte >25Gy ####
highdoses=' '
newline=1
for dummy in range(len(dose_matrix)):
    if dummy>3 and newline%6==0: 
        highdoses=highdoses+'\n '
        newline+=1
    if dose_matrix[dummy]>25: 
        highdoses=highdoses+str(np.int(np.round(otp_dose_values[dummy],0)))+'/'+str(np.int(np.round(dose_matrix[dummy],0)))+'Gy, '
        newline+=1
highdoses=highdoses[:-2]

cath_max_dose_5mm=cath_max_dose_5mm*corr_factor

#################################################
##### output ####################################
#################################################
err_diffs=[]
for dummy in range(len(check_points)):
    if otp_dose_values[dummy]<8:err_diffs.append(otp_dose_values[dummy]-dose_matrix[dummy])

otp_max=str(np.round(np.max(otp_dose_values),1))
pyth_max=str(np.round(np.max(dose_matrix),1))
otp_min=str(np.round(np.min(otp_dose_values),1))
pyth_min=str(np.round(np.min(dose_matrix),1))


gamma_5p=0
gamma_4p=0
gamma_3p=0
gamma_2p=0
gamma_1p=0
gamma_eval=0
for dummy in range(len(check_points)):
    if dose_matrix[dummy]<14: gamma_eval+=1
    if dose_matrix[dummy]<14 and np.abs(1-dose_matrix[dummy]/otp_dose_values[dummy])>0.05: gamma_5p+=1
    if dose_matrix[dummy]<14 and np.abs(1-dose_matrix[dummy]/otp_dose_values[dummy])>0.04: gamma_4p+=1
    if dose_matrix[dummy]<14 and np.abs(1-dose_matrix[dummy]/otp_dose_values[dummy])>0.03: gamma_3p+=1
    if dose_matrix[dummy]<14 and np.abs(1-dose_matrix[dummy]/otp_dose_values[dummy])>0.02: gamma_2p+=1
    if dose_matrix[dummy]<14 and np.abs(1-dose_matrix[dummy]/otp_dose_values[dummy])>0.01: gamma_1p+=1


######################################################################################
############# do plotting ############################################################

fig,axs=plt.subplots(3,2,sharey=False, tight_layout=True)
gs = axs[0, 0].get_gridspec()
for ax in axs[0:, 0]:
    ax.remove()
    axbig = fig.add_subplot(gs[0:, 0])
#    axbig.annotate('Big Axes \nGridSpec[0:, 0]', (0.1, 0.5),
#                   xycoords='axes fraction', va='center')
axbig.axis('off')

output_patient_data=' Pat: '+patname+' ['+patid+']'+'\n'+' Plan Datum: '+source_data3+' '+source_data4+'\n'

if len(highdoses)==0:
    output_auswertung=str(' Dosispunkte mit <5%-Abweichung: '+str(np.round(100*(1-gamma_5p/gamma_eval),1))+'%'+'\n'+
                          ' Dosispunkte mit <4%-Abweichung: '+str(np.round(100*(1-gamma_4p/gamma_eval),1))+'%'+'\n'+
                          ' Dosispunkte mit <3%-Abweichung: '+str(np.round(100*(1-gamma_3p/gamma_eval),1))+'%'+'\n'+
                          ' Dosispunkte mit <2%-Abweichung: '+str(np.round(100*(1-gamma_2p/gamma_eval),1))+'%'+'\n'+
                          ' Dosispunkte mit <1%-Abweichung: '+str(np.round(100*(1-gamma_1p/gamma_eval),1))+'%'+'\n\n'+
                          ' (*): -Berechnung auf Basis aller Punkte <14Gy: '+str(gamma_eval)+'/'+str(len(check_points))+'\n'+
                          '       -kein Gamma-Kriterium, Dosisvergleich OTP- vs. Python-'+'\n'+
                          '        Berechnung an exakt derselben geometrischen Position'+'\n\n'+
                          '\n\n Dosismaximum: '+otp_max+'Gy (OTP), '+pyth_max+'Gy (Python)'+'\n'+
                          ' Dosis größer 25Gy OTP/Python: keine')
else:    
    output_auswertung=str(' Dosispunkte mit <5%-Abweichung: '+str(np.round(100*(1-gamma_5p/gamma_eval),1))+'%'+'\n'+
                      ' Dosispunkte mit <4%-Abweichung: '+str(np.round(100*(1-gamma_4p/gamma_eval),1))+'%'+'\n'+
                      ' Dosispunkte mit <3%-Abweichung: '+str(np.round(100*(1-gamma_3p/gamma_eval),1))+'%'+'\n'+
                      ' Dosispunkte mit <2%-Abweichung: '+str(np.round(100*(1-gamma_2p/gamma_eval),1))+'%'+'\n'+
                      ' Dosispunkte mit <1%-Abweichung: '+str(np.round(100*(1-gamma_1p/gamma_eval),1))+'%'+'\n\n'+                  
                      ' (*): -Berechnung auf Basis aller Punkte <14Gy: '+str(gamma_eval)+'/'+str(len(check_points))+'\n'+
                      '       -kein Gamma-Kriterium, Dosisvergleich OTP- vs. Python-'+'\n'+
                      '        Berechnung an exakt derselben geometrischen Position'+'\n\n'+
                      '\n\n Dosismaximum: '+otp_max+'Gy (OTP), '+pyth_max+'Gy (Python)'+'\n'+
                      ' Dosis größer 25Gy OTP/Python: \n'+str(highdoses))

avg_rad_time=str(np.round(np.average(rad_time_positions),2))
output_kanalinfos=str(' Anzahl Kanäle: '+str(CN)+'\n'+' Kanal-Zeiten: '+text_cath_time+'\n'+
             ' ges. Strahlzeit: '+str(np.round(rad_time/60,1))+'min'+'\n'+
             ' Durchschnittliche Strahlzeit pro Position: '+avg_rad_time+'s')

infotext=str('[python-code v0.1 vom 11.3.2021, Medizinische Physik, Klinik für Strahlentherapie,'+'\n'+
             'Universitätsmedizin  Rostock, Südring 75, 18059 Rostock, verwendet anisotropie.txt'+'\n'+
             'aus OTP-Brachy vom 20.2.2021, Programm berechnet 3D Dosis neu anhand'+'\n'+ 
             'Quellenpositionen und Strahlzeiten aus RTPLAN-DICOM Datei und vergleicht'+'\n'+
             'diese mit den Dosispunkten aus der RTDOSE-DICOM Datei, die hier berechnete'+'\n'+
             'Dosis wurde skaliert auf die Berechnung aus OTP-Brachy, es erfolgt somit nur'+'\n'+
             'eine reine Kontrolle der 3D Dosisberechnung, nicht aber der tatsächlichen'+'\n'+
             'aktuellen Quellenaktivität]')

OUTPUT=str('\n'+
           '------------------------------------------------------------------------------------------------------------------'+'\n'+
           'Validierung der OTP-Brachy-3D-Dosisberechnung'+'\n'+
           '------------------------------------------------------------------------------------------------------------------'+'\n'+
           output_patient_data+'\n\n'
           '----------------------------------------------------------------------------------------------------'+'\n'+
           ' Auswertung (*):'+'\n'+
           '----------------------------------------------------------------------------------------------------'+'\n'+
           output_auswertung+'\n\n\n'+
           '----------------------------------------------------------------------------------------------------'+'\n'+
           ' Kanäle:'+'\n'+
           '----------------------------------------------------------------------------------------------------'+'\n'+
           output_kanalinfos)

axbig.text(0.05,1,OUTPUT,horizontalalignment='left',verticalalignment='top',fontsize=5.5)

axbig.text(0.1,0,infotext,fontsize=4.5)

####################################################################################
######################### HISTOGRAM PLOTS ##########################################
####################################################################################


axs[2,1].hist(err_diffs,bins=np.arange(-2, 2 + 0.01, 0.005))
axs[2,1].set_xlim([-0.15, 0.15])
axs[2,1].set_xlabel('$Dosis_{OTPBrachy} - Dosis_{Python} [Gy]$',fontsize=6)
axs[2,1].set_ylabel('$Dosispunkte (insg.: '+str(len(check_points))+')$',fontsize=6)
axs[2,1].set_title('Dosisdifferenz OTP-Python',fontsize=6)

if np.float(otp_max)>15: otp_max=15
bin_length=(np.float(otp_max)-np.float(otp_min))/30

axs[0,1].hist(otp_dose_values,bins=np.arange(0, 20 + 0.05, bin_length))
axs[0,1].set_xlim([np.float(otp_min)-0.2, np.float(otp_max)+0.2])
axs[0,1].set_xlabel('$Dosis_{OTPBrachy} [Gy]$',fontsize=6)
axs[0,1].set_ylabel('$Dosispunkte (insg.: '+str(len(check_points))+')$',fontsize=6)
axs[0,1].set_title('Dosis OTP-Berechnung',fontsize=6)

#axs[1,1].hist(dose_matrix,bins=np.arange(0, 20 + 0.05, bin_length))
axs[1,1].hist(dose_matrix_3d_all_points,bins=np.arange(0, 20 + 0.05, bin_length))
axs[1,1].set_xlim([np.float(otp_min)-0.2, np.float(otp_max)+0.2])
axs[1,1].set_xlabel('$Dosis_{Python} [Gy]$',fontsize=6)
axs[1,1].set_ylabel('$Dosispunkte (insg.: '+str(len(check_points))+')$',fontsize=6)
axs[1,1].set_title('Dosis Neuberechnung Python',fontsize=6)

plt.savefig('BrachyVerify_PatID'+patid+'_'+source_data3+'.pdf', bbox_inches='tight')



