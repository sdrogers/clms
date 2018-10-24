def MSplot(MStype,Data_mz_MS1=None,Data_intensity_MS1=None,Data_mz_MS2=None,Data_intensity_MS2=None,width=1,index=None,xlabel=0,ylabel=0,title=0,ylim=0, xlim=0,limits="Full",colours=None):
    import matplotlib.pyplot as plt
    import sys
    
    if MStype == "MS1":
        plt.bar(Data_mz_MS1, Data_intensity_MS1,width=width)
        if xlabel==0:
            plt.xlabel("m/z (MS1)")
        else:
            plt.xlabel(xlabel)
        
        if ylabel==0:
            plt.ylabel("Intensity")
        else:
            plt.ylabel(ylabel)
        
        if title==0:
            plt.title("MS1 Plot")
        else:
            plt.title(title)

        if xlim==0:
            plt.xlim([min(Data_mz_MS1)*0.95,max(Data_mz_MS1)*1.05])
        else:
            plt.xlim(xlim)
            
        if ylim==0:
            plt.ylim([min(Data_intensity_MS1)*0.95,max(Data_intensity_MS1)*1.05])
        else:
            plt.ylim(ylim)
            
        plt.show()        
    elif MStype == "MS2":
        if index is None :
            sys.exit("'index' must be defined for MS2")
            
        if colours is None:
            plt.bar(Data_mz_MS2[index], Data_intensity_MS2[index],width=width)
        else:
            p=plt.bar(Data_mz_MS2[index], Data_intensity_MS2[index],width=width)
            for i in range(0,len(Data_mz_MS2[index])):
                p[i].set_color(colours[i])
        
        if xlabel==0:
            plt.xlabel("m/z (MS2)")
        else:
            plt.xlabel(xlabel)
        
        if ylabel==0:
            plt.ylabel("Intensity")
        else:
            plt.ylabel(ylabel)
        
        if title==0:
            plt.title("MS2 Plot")
        else:
            plt.title(title)

        if xlim==0:
            plt.xlim([min(sum(Data_mz_MS2,[]))*0.95,max(sum(Data_mz_MS2,[]))*1.05])
        else:
            plt.xlim(xlim)
            
        if ylim==0:
            plt.ylim([min(sum(Data_intensity_MS2,[]))*0.95,max(sum(Data_intensity_MS2,[]))*1.05])
        else:
            plt.ylim(ylim)
            
        plt.show()
    else:
        sys.exit("Use either 'MS1' or 'MS2' for MStype")