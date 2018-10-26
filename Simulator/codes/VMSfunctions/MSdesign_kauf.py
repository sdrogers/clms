def MSdesign_kauf(Data,MS1_range_bottom,MS1_range_top,num_windows,window_type="Even"):
    import sys
    import numpy as np
    if (num_windows%2)!=0 or num_windows==0:
        sys.exit("'num_windows' must be even and non zero")
    else:
        if num_windows == 2 or num_windows == 4:
            print("This has not been implemented fully yet")
            #print("Design has not been optimised for this number of windows yet! Still works, but less effectively.")
        if num_windows > 8:
            print("This design may not be compatible with Mass Spec")
        if window_type=="Even": 
            scan_bottom_initial=[[MS1_range_bottom]]
            scan_top_initial=[[MS1_range_top]]
            scan_bottom_secondary=np.arange(MS1_range_bottom, MS1_range_top, (MS1_range_top-MS1_range_bottom)/num_windows).tolist()
            scan_bottom_secondary=[[x] for x in scan_bottom_secondary]
            scan_top_secondary=np.arange(MS1_range_bottom, MS1_range_top, (MS1_range_top-MS1_range_bottom)/num_windows).tolist()
            scan_top_secondary=[[x+(MS1_range_top-MS1_range_bottom)/num_windows] for x in scan_top_secondary]
            extra_scans_bottom = [[],[],[],[]]
            extra_scans_top = [[],[],[],[]]
            for i in range(0,int(num_windows),2):
                range_both = np.arange(sum(scan_bottom_secondary[i]),sum(scan_top_secondary[i+1]), (sum(scan_top_secondary[i+1])-sum(scan_bottom_secondary[i]))/16).tolist()
                range_bottom=[x for x in range_both]
                range_top = [x+(sum(scan_top_secondary[i+1])-sum(scan_bottom_secondary[i]))/16 for x in range_both] 
                # first extra scan
                extra_scans_bottom[0].extend([range_bottom[4]])
                extra_scans_top[0].extend([range_top[11]])
                # second extra scan
                extra_scans_bottom[1].extend([range_bottom[2],range_bottom[10]])
                extra_scans_top[1].extend([range_top[5],range_top[13]])
                # third extra scan
                extra_scans_bottom[2].extend([range_bottom[1],range_bottom[9]])
                extra_scans_top[2].extend([range_top[2],range_top[10]])
                # fouth extra scan
                extra_scans_bottom[3].extend([range_bottom[5],range_bottom[13]])
                extra_scans_top[3].extend([range_top[6],range_top[14]])
            scan_bottom_secondary.extend(extra_scans_bottom)
            scan_top_secondary.extend(extra_scans_top)
            MSdesign_result={'scan_bottom_initial':scan_bottom_initial,'scan_top_initial':scan_top_initial,'scan_bottom_secondary':scan_bottom_secondary,'scan_top_secondary':scan_top_secondary}
            return(MSdesign_result)
        elif window_type=="Percentile":
            scan_bottom_initial=[[MS1_range_bottom]]
            scan_top_initial=[[MS1_range_top]]
            scan_bottom_secondary = np.percentile(Data["mz_MS1"],np.arange(MS1_range_bottom, MS1_range_top, (MS1_range_top-MS1_range_bottom)/num_windows))
            scan_bottom_secondary=[[x] for x in scan_bottom_secondary]
            scan_bottom_secondary[0][0] = scan_bottom_secondary[0][0]*0.99 # adjustment for 0th percentile being min value
            x=np.arange(MS1_range_bottom, MS1_range_top+(MS1_range_top-MS1_range_bottom)/num_windows, (MS1_range_top-MS1_range_bottom)/num_windows)
            x1=x[1:len(x)]
            scan_top_secondary=np.percentile(Data["mz_MS1"],x1)
            scan_top_secondary=[[x] for x in scan_top_secondary]
            scan_top_secondary[-1][0]=scan_top_secondary[-1][0]*1.01 # adjustment for 100th percentile being max value
            extra_scans_bottom = [[],[],[],[]]
            extra_scans_top = [[],[],[],[]]
            for i in range(0,int(num_windows),2):
                scan_bottom_secondary2=np.arange(MS1_range_bottom, MS1_range_top, (MS1_range_top-MS1_range_bottom)/num_windows).tolist()
                scan_bottom_secondary2=[[x] for x in scan_bottom_secondary2]
                scan_top_secondary2=np.arange(MS1_range_bottom, MS1_range_top, (MS1_range_top-MS1_range_bottom)/num_windows).tolist()
                scan_top_secondary2=[[x+(MS1_range_top-MS1_range_bottom)/num_windows] for x in scan_top_secondary2]   
                range_both = np.arange(sum(scan_bottom_secondary2[i]),sum(scan_top_secondary2[i+1]), (sum(scan_top_secondary2[i+1])-sum(scan_bottom_secondary2[i]))/16).tolist()
                range_bottom=[x for x in range_both]
                range_top = [x+(sum(scan_top_secondary2[i+1])-sum(scan_bottom_secondary2[i]))/16 for x in range_both]
                range_top2=np.percentile(Data["mz_MS1"],range_top)
                range_bottom2=np.percentile(Data["mz_MS1"],range_bottom)
                if i==0: # I dont think these two if statements are necesary as we dont have use locations 0 and 15
                    range_bottom2[0]=range_bottom2[0]*0.99
                if i==(int(num_windows)-1):
                    range_top2[-1]=range_top2[-1]*1.01
                # first extra scan
                extra_scans_bottom[0].extend([range_bottom2[4]])
                extra_scans_top[0].extend([range_top2[11]])
                # second extra scan
                extra_scans_bottom[1].extend([range_bottom2[2],range_bottom2[10]])
                extra_scans_top[1].extend([range_top2[5],range_top2[13]])
                # third extra scan
                extra_scans_bottom[2].extend([range_bottom2[1],range_bottom2[9]])
                extra_scans_top[2].extend([range_top2[2],range_top2[10]])
                # fouth extra scan
                extra_scans_bottom[3].extend([range_bottom2[5],range_bottom2[13]])
                extra_scans_top[3].extend([range_top2[6],range_top2[14]])
            scan_bottom_secondary.extend(extra_scans_bottom)
            scan_top_secondary.extend(extra_scans_top)
            MSdesign_result={'scan_bottom_initial':scan_bottom_initial,'scan_top_initial':scan_top_initial,'scan_bottom_secondary':scan_bottom_secondary,'scan_top_secondary':scan_top_secondary}
            return(MSdesign_result)   
        else:
            sys.exit("window_type must be either 'Even' or 'Percentile' 1")
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            