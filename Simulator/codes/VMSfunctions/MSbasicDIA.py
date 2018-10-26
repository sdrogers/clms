def MSbasicDIA (Data,MS1_range_bottom,MS1_range_top,num_windows,window_type="Even"):
    import numpy as np
    from functions import MSmultiscan
    from functions import MSentropy
    # initial scan
    scan_bottom_initial=[[MS1_range_bottom]]
    scan_top_initial=[[MS1_range_top]]
    # secondary scans
    if window_type=="Even":
        scan_bottom_secondary=np.arange(MS1_range_bottom, MS1_range_top, (MS1_range_top-MS1_range_bottom)/num_windows).tolist()
        scan_bottom_secondary=[[x] for x in scan_bottom_secondary]
        scan_top_secondary=np.arange(MS1_range_bottom, MS1_range_top, (MS1_range_top-MS1_range_bottom)/num_windows).tolist()
        scan_top_secondary=[[x+(MS1_range_top-MS1_range_bottom)/num_windows] for x in scan_top_secondary]
    elif window_type=="Percentile":
        scan_bottom_secondary = np.percentile(Data.mz_ms1,np.arange(MS1_range_bottom, MS1_range_top, (MS1_range_top-MS1_range_bottom)/num_windows))
        scan_bottom_secondary=[[x] for x in scan_bottom_secondary]
        scan_bottom_secondary[0][0] = scan_bottom_secondary[0][0]*0.99 # adjustment for 0th percentile being min value
        x=np.arange(MS1_range_bottom, MS1_range_top+(MS1_range_top-MS1_range_bottom)/num_windows, (MS1_range_top-MS1_range_bottom)/num_windows)
        x1=x[1:len(x)]
        scan_top_secondary=np.percentile(Data.mz_ms1,x1)
        scan_top_secondary=[[x] for x in scan_top_secondary]
        scan_top_secondary[-1][0]=scan_top_secondary[-1][0]*1.01 # adjustment for 100th percentile being max value
    else:
        sys.exit("window_type must be 'Even' or 'Percentile'")    
    # do the scans
    DIA_result = Data.MS_MultiScan(scan_bottom_initial,scan_top_initial,scan_bottom_secondary,scan_top_secondary)
    # calculate ranges
    DIA_result2 = DIA_result["Scan Results"][0]
    DIA_lower = []
    DIA_upper = []
    for i in range(0,len(DIA_result2)):
        j=1
        found=False
        while found==False:
            if DIA_result["Scan Results"][j+1][i]==1:
                found=True
                DIA_lower.append(DIA_result["Scan Locations"][j][0])
                DIA_upper.append(DIA_result["Scan Locations"][j][1])
            else:
                j+=1
    # calculate entropy
    entropy=MSentropy.MSentropy(Data,DIA_lower,DIA_upper)
    DIA_results = {'MS2 Locations':DIA_result2,'DIA_lower':DIA_lower, 'DIA_upper':DIA_upper, 'entropy':entropy, 'Scan Results':DIA_result["Scan Results"],'Scan Locations':DIA_result["Scan Locations"]}
    return(DIA_results)
    
    
    