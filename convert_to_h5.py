#!/usr/bin/env python3

import sys
import getopt
import h5py
import ROOT
import numpy as np
import os
from tqdm import tqdm

ROOT.gSystem.Load(f"{os.environ['WCSIM_BUILD_DIR']}/lib/libWCSimRoot.so")

def usage():
    '''Demo script to convert wcsim.root to h5
    '''
    print ("Demo script to convert wcsim.root to h5")
    print ("Usage:")
    print ("convert_to_h5.py [-h] [-f <wcsim_file_to_convert>] [-q <fiTQun_file>] [-o <output_filename>] [-m <mask_file>]")
    print ("")
    print ("Options:")
    print ("-h, --help: prints help message")
    print ("-f, --file=<file_to_convert>")
    print ("-q, --fqfile=<fiTQun_file>")
    print ("-o, --out=<output_filename>")
    print ("")

def convert_to_h5():

    fname = None
    fqname= None
    mask_file_name = "maskFile_1766.txt"
    fout = 'out.h5'

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:o:q:m:",
                                   ["help", "file=" , "fqfile=", "out=", "mask="])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
        sys.exit(2)

    for opt, val in opts:
        if (opt in ("-h", "--help")):
            usage()
            sys.exit()
        if (opt in ("-f", "--file")):
            fname = val.strip()
            print(f"fname =  {fname}")
        if (opt in ("-q", "--fqfile")):
            fqname = val.strip()
            print(f"fqname =  {fqname}")
        if (opt in ("-o", "--out")):
            fout = val.strip()
        if (opt in ("-m", "--mask")):
            mask_file_name = val.strip()
            print(f"mask_file_name =  {mask_file_name}")

    if fname == None:
        print("Missing wcsim input file!!")
        usage()
        sys.exit(2)

    if fqname == None:
        print("Missing fiTQun file!!")
        usage()
        sys.exit(2)

    chain = ROOT.TChain("wcsimT")
    chain.Add(fname)

    max_r = 3075.926/2.0
    max_z = 2714.235/2.0

    single_file_name = chain.GetFile().GetName()
    f = ROOT.TFile.Open(single_file_name)
    if not f or not f.IsOpen():
        print(f"Error, could not open input file: {single_file_name}")
        return -1
    geotree = f.Get("wcsimGeoT")
    print(f"Geotree has {geotree.GetEntries()} entries")
    if geotree.GetEntries() == 0:
        exit(9)
    geotree.GetEntry(0)
    geo = geotree.wcsimrootgeom
    npmts = geo.GetWCNumPMT()
    print(f"geo has {npmts} PMTs")
    pmt_CylLoc = [0]*npmts
    for i in range(npmts):
        pmt = geo.GetPMT(i)
        pmt_CylLoc[i] = pmt.GetCylLoc()

    # Friend chain (fiTQun)
    chainFQ = ROOT.TChain("fiTQun")
    chainFQ.Add(fqname)
    chain.AddFriend(chainFQ)

    # PMT bad channel list
    mask_list = [False] * npmts
    with open(mask_file_name, 'r') as maskFile:
        for line in maskFile:
            cab_id = int(line.strip())
            mask_list[cab_id] = True
            # print(f"cab_id {cab_id} is masked")

    nevents = chain.GetEntries()
    print("number of entries in the tree: " + str(nevents))
    chain.GetEvent(0)
    event = chain.wcsimrootevent

    vertex = np.zeros((nevents,3),dtype=np.float32)
    direction = np.zeros((nevents,3),dtype=np.float32)
    entrance_pos = np.zeros((nevents,3),dtype=np.float32)
    exit_pos = np.zeros((nevents,3),dtype=np.float32)
    is_top_down = np.ones(nevents,dtype=bool)
    momentum = np.zeros(nevents,dtype=np.float32)
    pmtQ = np.zeros((nevents,npmts),dtype=np.float32)
    pmtT = np.zeros((nevents,npmts),dtype=np.float32)
    fq_pos = np.zeros((nevents,4),dtype=np.float32)
    fq_entrance_pos = np.zeros((nevents,3),dtype=np.float32)
    fq_exit_pos = np.zeros((nevents,3),dtype=np.float32)
    fq_direction = np.zeros((nevents,3),dtype=np.float32)
    fq_momentum = np.zeros(nevents,dtype=np.float32)
    pass_selection = np.ones(nevents,dtype=bool)

    for i in tqdm(range(nevents)):
        event.ReInitialize()
        chain.GetEntry(i)
        event = chain.wcsimrootevent
        trigger = event.GetTrigger(0)

        tracks = trigger.GetTracks()
        primary_track = [t for t in tracks if (t.GetId() == 1)]
        vertex[i] = np.array([primary_track[0].GetStart(0)*10, -primary_track[0].GetStart(2)*10, primary_track[0].GetStart(1)*10], dtype=np.float32) # in mm
        direction[i] = np.array([primary_track[0].GetDir(0), -primary_track[0].GetDir(2), primary_track[0].GetDir(1)], dtype=np.float32)
        momentum[i] = primary_track[0].GetP()

        bp = primary_track[0].GetBoundaryPoints()
        bt = primary_track[0].GetBoundaryTypes()
        if len(bp)>=4:
            entrance_pos[i] = np.array([bp[0][0], -bp[0][2], bp[0][1]], dtype=np.float32) # in mm
            if bt[0]==1 and entrance_pos[i][2]<max_z:  # enter through blacksheet
                is_top_down[i] = False
            elif bt[0]==2 and entrance_pos[i][2]<max_z-100:  # enter through mPMT
                is_top_down[i] = False
            last = len(bp)-1
            exit_pos[i] = np.array([bp[last][0], -bp[last][2], bp[last][1]], dtype=np.float32) # in mm
            if bt[last]==1 and exit_pos[i][2]>-max_z:  # exit through blacksheet
                is_top_down[i] = False
            elif bt[last]==2 and exit_pos[i][2]>-max_z+100:  # exit through mPMT
                is_top_down[i] = False

        # reject short track
        if (np.linalg.norm(entrance_pos[i] - exit_pos[i])<500):
            pass_selection[i] = False
            is_top_down[i] = False

        nhits = 0
        totQ = 0.0
        pmt_hit = [0.0, 0.0, 0.0]
        for hit in trigger.GetCherenkovDigiHits():
            if hit.GetT()<100:
                pmt_id = hit.GetTubeId() - 1
                if mask_list[pmt_id]:
                    continue
                q = hit.GetQ()
                pmtQ[i][pmt_id] = pmtQ[i][pmt_id]+q
                if pmtT[i][pmt_id]<=0:
                    pmtT[i][pmt_id] = hit.GetT()
                cycloc = pmt_CylLoc[pmt_id]
                pmt_hit[cycloc] += q
                nhits += 1
                totQ += q
        
        fq1rpos = chain.fq1rpos  # expected nested array: [eventIndex][fitIndex][xyz]
        fq1rdir = chain.fq1rdir
        fq1rmom = chain.fq1rmom
        fq1rt0  = chain.fq1rt0

        fqtopdown = True
        fqVtx = ROOT.TVector3(fq1rpos[2*3+0]*10, -fq1rpos[2*3+2]*10, fq1rpos[2*3+1]*10)
        fqDir = ROOT.TVector3(fq1rdir[2*3+0],    -fq1rdir[2*3+2],    fq1rdir[2*3+1])

        if fqDir.Z() == 0:
            fqtopdown = False
        else:
            topcapdist = (max_z - fqVtx.Z()) / fqDir.Z()
            topPt = fqVtx + topcapdist * fqDir
            if topPt.Perp() > max_r:
                # move along direction until hits barrel or exit z-range
                while topPt.Perp() > max_r and abs(topPt.Z()) <= max_z:
                    topPt = topPt + fqDir
                fqtopdown = False

            botcapdist = (-max_z - fqVtx.Z()) / fqDir.Z()
            botPt = fqVtx + botcapdist * fqDir
            if botPt.Perp() > max_r:
                while botPt.Perp() > max_r and abs(botPt.Z()) <= max_z:
                    botPt = botPt - fqDir
                fqtopdown = False

        fq_pos[i] = np.array([fqVtx.X(),fqVtx.Y(),fqVtx.Z(),fq1rt0[2]], dtype=np.float32)
        fq_entrance_pos[i] = np.array([topPt.X(),topPt.Y(),topPt.Z()], dtype=np.float32)
        fq_exit_pos[i] = np.array([botPt.X(),botPt.Y(),botPt.Z()], dtype=np.float32)
        fq_direction[i] = np.array([fqDir.X(),fqDir.Y(),fqDir.Z()], dtype=np.float32)
        fq_momentum[i] = fq1rmom[2]

        if totQ == 0:
            pass_selection[i] = False
        else:
            if pmt_hit[0]/totQ > 0.07:
                pass_selection[i] = False
            if pmt_hit[1]/totQ < 0.45 or pmt_hit[1]/totQ > 0.7:
                pass_selection[i] = False
            if pmt_hit[2]/totQ < 0.25 or pmt_hit[2]/totQ > 0.5:
                pass_selection[i] = False
            if nhits < 700:
                pass_selection[i] = False
            if fqDir.Z() > 0:
                pass_selection[i] = False
            if not fqtopdown:
                pass_selection[i] = False

    print('vertex.shape = ',vertex.shape)
    print('direction.shape = ',direction.shape)
    print('momentum.shape = ',momentum.shape)
    print('pmtQ.shape = ',pmtQ.shape)
    print('# selected = ',np.sum(pass_selection))
    print('# true top-down = ',np.sum(is_top_down))
    print('# (selected and true top-down)  = ',np.sum(pass_selection*is_top_down))


    with h5py.File(fout,mode='w') as h5fw:
        h5fw.create_dataset('vertex', data=vertex)
        h5fw.create_dataset('direction', data=direction)
        h5fw.create_dataset('entrance_pos', data=entrance_pos)
        h5fw.create_dataset('exit_pos', data=exit_pos)
        h5fw.create_dataset('is_top_down', data=is_top_down)
        h5fw.create_dataset('momentum', data=momentum)
        h5fw.create_dataset('pmtQ', data=pmtQ)
        h5fw.create_dataset('pmtT', data=pmtT)
        h5fw.create_dataset('fq_pos', data=fq_pos)
        h5fw.create_dataset('fq_entrance_pos', data=fq_entrance_pos)
        h5fw.create_dataset('fq_exit_pos', data=fq_exit_pos)
        h5fw.create_dataset('fq_direction', data=fq_direction)
        h5fw.create_dataset('fq_momentum', data=fq_momentum)
        h5fw.create_dataset('pass_selection', data=pass_selection)

if __name__ == '__main__':
    convert_to_h5()