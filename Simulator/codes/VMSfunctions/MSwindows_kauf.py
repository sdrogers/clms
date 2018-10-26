def MSwindows_kauf (Data,DIA_result,window_type="Even"):
    from functions import MSentropy
    import sys
    import numpy as np
    import math
    if window_type=="Even":
        template1= [[0,0,0,0],[0,0,1,0],[0,1,1,0],[0,1,0,0],[1,1,0,0],[1,1,0,1],[1,0,0,1],[1,0,0,0]]
        template2=[[1,0,0,0],[1,0,1,0],[1,1,1,0],[1,1,0,0],[0,1,0,0],[0,1,0,1],[0,0,0,1],[0,0,0,0]]
        DIA_result2 = DIA_result["Scan Results"][0]
        DIA_lower = []
        DIA_upper = []
        for i in range(0,len(DIA_result2)):
            j=1
            found=False
            while found==False:
                if DIA_result["Scan Results"][j+1][i]==1:
                    found=True
                    initial_lower=DIA_result["Scan Locations"][j][0]
                    initial_higher=DIA_result["Scan Locations"][j][1]
                    final4extract=[DIA_result["Scan Results"][10][i],DIA_result["Scan Results"][11][i],DIA_result["Scan Results"][12][i],DIA_result["Scan Results"][13][i]]
                    range_both = np.arange(initial_lower,initial_higher, (initial_higher-initial_lower)/8).tolist()
                    range_bottom=[x for x in range_both]
                    range_top = [x+(initial_higher-initial_lower)/8 for x in range_both] 
                    check =False
                    c=0
                    if math.floor(j/2)==(j/2):
                        while check ==False:
                            if template2[c]==final4extract:
                                check = True
                                DIA_lower.append(range_bottom[c])
                                DIA_upper.append(range_top[c])
                            c+=1
                    else:    
                        while check ==False:
                            if template1[c]==final4extract:
                                check = True
                                DIA_lower.append(range_bottom[c])
                                DIA_upper.append(range_top[c])
                            c+=1
                else:
                    j+=1
        entropy=MSentropy.MSentropy(Data,DIA_lower,DIA_upper)
        DIA_results = {'MS2 Locations':DIA_result2,'DIA_lower':DIA_lower, 'DIA_upper':DIA_upper, 'entropy':entropy, 'Scan Results':DIA_result["Scan Results"],'Scan Locations':DIA_result["Scan Locations"]}
        return(DIA_results)
    elif window_type=="Percentile":
        template1= [[0,0,0,0],[0,0,1,0],[0,1,1,0],[0,1,0,0],[1,1,0,0],[1,1,0,1],[1,0,0,1],[1,0,0,0]]
        template2=[[1,0,0,0],[1,0,1,0],[1,1,1,0],[1,1,0,0],[0,1,0,0],[0,1,0,1],[0,0,0,1],[0,0,0,0]]
        DIA_result2 = DIA_result["Scan Results"][0]
        DIA_lower = []
        DIA_upper = []
        for i in range(0,len(DIA_result2)):
            j=1
            found=False
            while found==False:
                if DIA_result["Scan Results"][j+1][i]==1:
                    found=True
                    initial_lower=DIA_result["Scan Locations"][j][0]
                    initial_higher=DIA_result["Scan Locations"][j][1]
                    final4extract=[DIA_result["Scan Results"][10][i],DIA_result["Scan Results"][11][i],DIA_result["Scan Results"][12][i],DIA_result["Scan Results"][13][i]]
                    bounds=np.percentile(Data["mz_MS1"],np.arange(0,100+100/64,100/64))
                    bounds[0]=bounds[0]*0.99
                    bounds[-1]=bounds[-1]*1.01
                    range_bottom=bounds[np.arange(int(np.where(bounds==initial_lower)[0]),int(np.where(bounds==initial_lower)[0])+8,1)]
                    range_top=bounds[np.arange(int(np.where(bounds==initial_lower)[0])+1,int(np.where(bounds==initial_lower)[0])+8+1,1)]
                    check =False
                    c=0
                    if math.floor(j/2)==(j/2):
                        while check ==False:
                            if template2[c]==final4extract:
                                check = True
                                DIA_lower.append(range_bottom[c])
                                DIA_upper.append(range_top[c])
                            c+=1
                    else:    
                        while check ==False:
                            if template1[c]==final4extract:
                                check = True
                                DIA_lower.append(range_bottom[c])
                                DIA_upper.append(range_top[c])
                            c+=1
                else:
                    j+=1
        entropy=MSentropy.MSentropy(Data,DIA_lower,DIA_upper)
        DIA_results = {'MS2 Locations':DIA_result2,'DIA_lower':DIA_lower, 'DIA_upper':DIA_upper, 'entropy':entropy, 'Scan Results':DIA_result["Scan Results"],'Scan Locations':DIA_result["Scan Locations"]}
        return(DIA_results)
    else:
        sys.exit("window_type must be either 'Even' or 'Percentile' 2")