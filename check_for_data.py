#!/usr/bin/python
import psutil
import pydicom as dcm
import os.path
import glob,os
import sys
import time
import subprocess
from os import path

importpath="/mnt/w/Verifikation_RS-Monaco/*"
backuppath="/mnt/w/Sicherung/BrachyVerify/"


run_check=glob.glob('currently_running.sys')

if len(run_check)==0:
    os.popen('touch currently_running.sys')

    patients_found_RP=[]
    patients_found_RD=[]
    patients_found_folder=[]
    patfolder=[]

    folders=glob.glob(importpath)
    for dummy in range(len(folders)):
        ftest_all_in_folder=glob.glob(folders[dummy]+'/*dcm')
        ftest_P=''
        ftest_D=''
        for file_runner in range(len(ftest_all_in_folder)):
            fcheck=str(ftest_all_in_folder[file_runner])
            fcheck=fcheck.split('/')
            fcheck=fcheck[-1]
            if fcheck.find('Monaco')==-1 and fcheck.find('RayStation')==-1:
                dcm_test=dcm.read_file(ftest_all_in_folder[file_runner],force=True)
                data_str=str(dcm_test)
                if data_str.find('BRACHY')>0 and dcm_test.Modality=='RTPLAN':
                   ftest_P=ftest_all_in_folder[file_runner]
                   print('found Pat: ',ftest_P)
                if data_str.find('MDS NORDION CALCULATION')>0 and dcm_test.Modality=='RTDOSE':
                   ftest_D=ftest_all_in_folder[file_runner]
                   print('found Pat: ',ftest_D)

        if len(ftest_P)>0 and len(ftest_D)>0:
            patients_found_RP.append(ftest_P)
            patients_found_RD.append(ftest_D)
            patients_found_folder.append(folders[dummy])

    if len(patients_found_folder)>0: print(patients_found_folder)

    for dummy in range(len(patients_found_folder)):
        patfolder_=patients_found_folder[dummy].split('/')
        patfolder.append(patfolder_[-1])
        mainpy=subprocess.Popen(['./BRACHY_VERIFY.py',str(patients_found_RP[dummy]),str(patients_found_RD[dummy])])
        mainpy.wait()

        renamer=[]
        renamer.append('')
        for rr in range(50):
            renamer.append('_'+str(rr+1))
        found_counter=0
        found_folder=1
        while found_folder==1:
            if os.path.isdir(backuppath+patfolder[dummy]+renamer[found_counter])==True:
                found_counter+=1
            else:
                found_folder=0

        os.popen('mv '+patients_found_folder[dummy]+' '+backuppath+'/'+patfolder[dummy]+renamer[found_counter])
        time.sleep(5)
        os.popen('cp *.pdf '+backuppath+patfolder[dummy]+renamer[found_counter]+'/')
        time.sleep(5)
        os.popen('mv *.pdf '+backuppath+'/')

    run_check=glob.glob('currently_running.sys')
    if len(run_check)==1: os.popen('rm currently_running.sys')
