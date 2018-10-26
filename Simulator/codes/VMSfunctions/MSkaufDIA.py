def MSkaufDIA(Data,MS1_range_bottom,MS1_range_top,num_windows,window_type="Even"):
    import sys
    import numpy as np
    from functions import MSmultiscan
    from functions import MSdesign_kauf
    from functions import MSwindows_kauf
    if window_type=="Even":          
        MSdesign_result=MSdesign_kauf.MSdesign_kauf(Data,MS1_range_bottom,MS1_range_top,num_windows,window_type=window_type)
    elif window_type=="Percentile":
        MSdesign_result=MSdesign_kauf.MSdesign_kauf(Data,MS1_range_bottom,MS1_range_top,num_windows,window_type=window_type)
    else:
        MSdesign_result=sys.exit("window_type must be either 'Even' or 'Percentile' 3")
    DIA_result = MSmultiscan.MSmultiscan(Data,scan_bottom_initial=MSdesign_result["scan_bottom_initial"],scan_top_initial=MSdesign_result["scan_top_initial"],scan_bottom_secondary=MSdesign_result["scan_bottom_secondary"],scan_top_secondary=MSdesign_result["scan_top_secondary"])
    if window_type=="Even":
        DIA_results=MSwindows_kauf.MSwindows_kauf(Data,DIA_result,window_type=window_type)
    elif window_type=="Percentile":
        DIA_results=MSwindows_kauf.MSwindows_kauf(Data,DIA_result,window_type=window_type)
    else:
        sys.exit("window_type must be either 'Even' or 'Percentile' 4")
    return(DIA_results)
    
    
    
    
    
    
    
    
    
    
    