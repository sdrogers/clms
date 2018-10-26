def MSentropy (Data,DIA_lower,DIA_upper):
    import math
    entropy_vec = []
    MS1=Data["mz_MS1"]
    for i in range(0,len(DIA_lower)):
        entropy_vec.extend([0])
        for j in range(0,len(MS1)):
            if DIA_lower[i] < MS1[j] and DIA_upper[i] > MS1[j]:
                entropy_vec[i] +=1
    entropy = sum([math.log(y) for y in entropy_vec])
    return(entropy)