import numpy as np
class MS_Data(object):
  def __init__(self,dataset,rt):
    self.mz_ms1 = []
    self.intensity_ms1 = []
    self.mz_ms2 = []
    self.intensity_ms2 = []
    self.ms1_object = []
    for i in range(0,len(dataset)):
        if(dataset[i].rt==rt):
            if(dataset[i].ms_level==1):
                if(dataset[i] in self.ms1_object):
                    self.mz_ms1[self.ms1_object.index(dataset[i])] = dataset[i].mz
                    self.intensity_ms1[self.ms1_object.index(dataset[i])] = dataset[i].intensity
                else:
                    self.mz_ms1.append(dataset[i].mz)
                    self.intensity_ms1.append(dataset[i].intensity)
                    self.mz_ms2.append(None)
                    self.intensity_ms2.append(None)
                    self.ms1_object.append(dataset[i]) 
            else:
                if(dataset[i].parent in self.ms1_object):
                    self.mz_ms2[self.ms1_object.index(dataset[i].parent)] = dataset[i].mz
                    self.intensity_ms2[self.ms1_object.index(dataset[i].parent)] = dataset[i].intensity
                else:
                    self.mz_ms1.append(None)
                    self.intensity_ms1.append(None)
                    self.mz_ms2.append(dataset[i].mz)
                    self.intensity_ms2.append(dataset[i].intensity)
                    self.ms1_object.append(dataset[i].parent) 
                    
  def MS_Plot(self,MS_type,Index=None,width=1,xlabel=0,ylabel=0,title=0,ylim=0,xlim=0,colours=None):
    MSplot(MS_type,self.mz_ms1,self.intensity_ms1,self.mz_ms2,self.intensity_ms2,index=Index,width=1,xlabel=0,ylabel=0,title=0,ylim=0,xlim=0,colours=None)   

  def MS_Scan(self,range_min,range_max,plot=False,colour=False):
    which_in_range = np.zeros(len(self.mz_ms1))
    for i in range(0,len(range_min)):
        which_in_range = np.add(which_in_range,np.concatenate(((np.array(self.mz_ms1)<=range_max[i]) & (np.array(self.mz_ms1)>range_min[i]))*1))
    which_in_range2 = np.array(which_in_range)==1
    mylist=np.array(range(0,len(self.mz_ms1)))[which_in_range2]
    mz_ms2_list = [self.mz_ms2[i] for i in mylist]
    intensity_ms2_list = [self.intensity_ms2[i] for i in mylist]
    return([mz_ms2_list,intensity_ms2_list])
    # NEED TO RETURN AN OBJECT HERE SOMEHOW???


    #if plot==True:
    #    intensity_ms2_list = [x for _,x in sorted(zip(mz_ms2_list,intensity_ms2_list))]
    #    mz_ms2_list.sort()
    #    if colour==True:
    #        available_colours=["b","g","r","c","m","y","k"]
    #        colours=[]
    #        for i in range(0,len(mz_ms2_list)):
    #            colours.append([available_colours[i%7]]*len(mz_ms2_list[i]))
    #        colours = sum(colours,[])
    #        MSplot(MStype="MS2",Data_mz_MS2=[sum(mz_ms2_list,[])],Data_intensity_ms2=[sum(intensity_ms2_list,[])],index=0,colours=colours)
    #    else:
    #        MSplot(MStype="MS2",Data_mz_MS2=[sum(mz_ms2_list,[])],Data_intensity_ms2=[sum(intensity_ms2_list,[])],index=0)

  def MS_MultiScan(self,scan_bottom_initial,scan_top_initial,scan_bottom_secondary,scan_top_secondary):
    scans = []
    scan_location = []
    scan_result = [[] for e in range(0,len(scan_bottom_initial))] 
    for i in range(0,len(scan_bottom_initial)):
        which_in_range = np.zeros(len(self.mz_ms1))
        for j in range(0,len(scan_bottom_initial[i])):
            which_in_range = np.add(which_in_range,np.concatenate(((np.array(self.mz_ms1)<=scan_top_initial[i][j]) & (np.array(self.mz_ms1)>scan_bottom_initial[i][j]))*1))
        which_in_range2 = np.array(which_in_range)==1
        mylist=np.array(range(0,len(self.mz_ms1)))[which_in_range2]
        mz_ms2_list = [self.mz_ms2[i] for i in mylist]
        scans.extend(sum(mz_ms2_list,[]))
        scan_location.append(sum([scan_bottom_initial[i],scan_top_initial[i]],[]))
        for j in range(0,len(scan_bottom_initial)):
            if j==i:
                scan_result[j].extend([1]*len(sum(mz_ms2_list,[])))
            else:
                scan_result[j].extend([0]*len(sum(mz_ms2_list,[])))
    if len(set(scans))!=len(scans):
        sys.exit("An area has been scanned twice in the initial set up phase")

    for i in range(0,len(scan_bottom_secondary)):
        which_in_range = np.zeros(len(self.mz_ms1))
        for j in range(0,len(scan_bottom_secondary[i])):
            which_in_range = np.add(which_in_range,np.concatenate(((np.array(self.mz_ms1)<=scan_top_secondary[i][j]) & (np.array(self.mz_ms1)>scan_bottom_secondary[i][j]))*1))
        which_in_range2 = np.array(which_in_range)==1
        mylist=np.array(range(0,len(self.mz_ms1)))[which_in_range2]
        mz_ms2_list = [self.mz_ms2[i] for i in mylist]
        if len(scan_bottom_secondary[i])==1:
            scan_location.append(sum([scan_bottom_secondary[i],scan_top_secondary[i]],[]))
        else:
            scan_location.append([scan_bottom_secondary[i],scan_top_secondary[i]])
        ms2_s = sum(mz_ms2_list,[])
        s_scan = [None]*len(scans)
        for j in range(0,len(scans)):
            if scans[j] in sum(mz_ms2_list,[]):
                s_scan[j]=1
            else:
                s_scan[j]=0
        scan_result.extend([s_scan])
    scans=[scans]
    scans.extend(scan_result)
    return({'Scan Results':scans, 'Scan Locations':scan_location})