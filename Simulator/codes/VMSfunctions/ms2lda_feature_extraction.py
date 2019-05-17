# coding=utf8
# Simon's attempts to make a single feature selection pipeline
from __future__ import print_function
# from Queue import PriorityQueue
import numpy as np
import sys, os
import re
import json
# sys.path.append('/Users/simon/git/efcompute')
# from ef_assigner import ef_assigner
# from formula import Formula
# from ef_constants import ATOM_MASSES, PROTON_MASS, ATOM_NAME_LIST

# Restructuring in December 2016
# Feature selection has three steps:
#  1. Loading a bunch of spectra into some standard format
#  2. Turning them into fragment and loss features
#  3. Making the corpus object


PROTON_MASS = 1.00727645199076
class MS1(object):
    def __init__(self,id,mz,rt,intensity,file_name,scan_number = None,single_charge_precursor_mass = None):
        self.id = id
        self.mz = mz
        self.rt = rt
        self.intensity = intensity
        self.file_name = file_name
        self.scan_number = scan_number
        if single_charge_precursor_mass:
            self.single_charge_precursor_mass = single_charge_precursor_mass
        else:
            self.single_charge_precursor_mass = self.mz
        self.name = "{}_{}".format(self.mz,self.rt)

    def __str__(self):
        return self.name






# ******************************
# ******************************
# ******************************
# LOADERS
# ******************************
# ******************************
# ******************************




# Abstract loader class
## Refactored Sep 21, 2017
## class LoadMZML, LoadMSP, LoadMGF will inhereit from class Loader
## Three sub-classes will implement their own *load_spectra* function based on different input ms2 files(mzml, msp, mgf)

## *load_spectra* functions are too long, refactor and split when having time
class Loader(object):
    def __init__(self,min_ms1_intensity = 0.0,peaklist = None,isolation_window = 0.5,mz_tol = 5,rt_tol=5.0,duplicate_filter_mz_tol = 0.5,duplicate_filter_rt_tol = 16,duplicate_filter = False,repeated_precursor_match = None,
                    min_ms1_rt = 0.0, max_ms1_rt = 1e6, min_ms2_intensity = 0.0,has_scan_id = False, rt_units = 'seconds',mz_col_name = 'mz', rt_col_name = 'rt', csv_id_col = None, id_field = None,name_field = None):
        self.min_ms1_intensity = min_ms1_intensity
        self.peaklist = peaklist
        self.isolation_window = isolation_window
        self.mz_tol = mz_tol
        self.rt_tol = rt_tol
        self.duplicate_filter = duplicate_filter
        self.duplicate_filter_mz_tol = duplicate_filter_mz_tol
        self.duplicate_filter_rt_tol = duplicate_filter_rt_tol
        self.min_ms1_rt = min_ms1_rt
        self.max_ms1_rt = max_ms1_rt
        self.min_ms2_intensity = min_ms2_intensity
        if repeated_precursor_match:
            self.repeated_precursor_match = repeated_precursor_match
        else:
            self.repeated_precursor_match = 2*self.isolation_window


        self.mz_col_name = mz_col_name
        self.rt_col_name = rt_col_name
        self.csv_id_col = csv_id_col
        self.rt_units = rt_units
        self.csv_id_col = csv_id_col
        self.id_field = id_field

        self.name_field = name_field # only works for msp - fix for metlin people

        if not self.mz_col_name:
            self.mz_col_name = 'mz'

    def __str__(self):
        return "loader class"

    def load_spectra(self,input_set):
        raise NotImplementedError("load spectra method must be implemented")


    # compute the parent masses
    # single_charge version is used for loss computation
    def _ion_masses(self,precursormass,int_charge):
        mul = abs(int_charge)
        parent_mass = precursormass*mul
        parent_mass -= int_charge*PROTON_MASS
        single_charge_precursor_mass = precursormass*mul
        if int_charge > 0:
            single_charge_precursor_mass -= (int_charge-1)*PROTON_MASS
        elif int_charge < 0:
            single_charge_precursor_mass += (mul-1)*PROTON_MASS
        else:
            # charge = zero - leave them ll the same
            parent_mass = precursormass
            single_charge_precursor_mass = precursormass
        return parent_mass,single_charge_precursor_mass


    ## method to interpret the ever variable charge
    ## field in the different formats
    ## should never fail now
    def _interpret_charge(self,charge):
        if not charge: # if it is none
            return 1
        try:
            if not type(charge) == str:
                charge = str(charge)

            ## add the meat here
            ## try removing any + signs
            charge = charge.replace("+", "")

            ## remove trailing minus signs
            if charge.endswith('-'):
                charge = charge[:-1]
                # move the minus to the front if it 
                # isn't already there
                if not charge.startswith('-'):
                    charge = '-' + charge
            ## turn into an int
            int_charge = int(charge)
            return int_charge
        except:
            int_charge = 1
        return int_charge
    ## modify peaklist function
    ## try to detect "featureid", store it in ms1_peaks used for in for mgf ms1 analysis
    ## ms1_peaks: [featid, mz,rt,intensity], featid will be None if "FeatureId" not exist
    def _load_peak_list(self):
        self.ms1_peaks = []
        self.user_cols_names = []
        with open(self.peaklist,'rU') as f:


            heads = f.readline()

            ## add this in case peaklist file is separated by ';'
            self.separator = ','
            if ';' in heads:
                self.separator = ';'

            tokens = heads.strip().split(self.separator)
            index = -1
            featid_index = None
            mz_col = None
            rt_col = None
            for i in range(len(tokens)):
                if tokens[i].lower() == self.mz_col_name.lower():
                    index = i
                elif self.csv_id_col and tokens[i].lower() == self.csv_id_col.lower():
                    featid_index = i
                # if tokens[i].lower() == "scans":
                #     featid_index = i
                if tokens[i].lower() in ['mass', 'mz']: # backwards compatibility
                    index = i
                #     break
                self.user_cols_names.append(tokens[i])

            ## if any sample names missing, use "Sample_*" to replace
            empty_sample_name_id = 0
            for i in range(index+2, len(tokens)):
                if not tokens[i]:
                    tokens[i] = "Sample_" + str(empty_sample_name_id)
                    empty_sample_name_id += 1

            self.sample_names = tokens[index+2:]

            for line in f:
                tokens_tuple= line.strip().split(self.separator, index+2)
                featid = None
                if featid_index != None:
                    featid = tokens_tuple[featid_index]
                mz = tokens_tuple[index]
                rt = float(tokens_tuple[index+1])
                if self.rt_units == 'minutes':
                    rt *= 60.0
                samples = tokens_tuple[index+2]
                # store (featid, mz,rt,intensity)

                ## record user defined index columns before "mass" column in peaklist file
                try:
                    self.ms1_peaks.append((featid, float(mz), float(rt), samples, tokens_tuple[:index]))
                except:
                    print("Failed on line: ")
                    print(line)

        # sort them by mass
        self.ms1_peaks = sorted(self.ms1_peaks,key = lambda x: x[1])
        print("Loaded {} ms1 peaks from {}".format(len(self.ms1_peaks),self.peaklist))

    ## read in peaklist file (.csv)
    ## ("..., mass, RT, samplename_1, samplename_2,..."), delimiter: '.
    ## find the most suitable ms1 hit
    ## then update ms1, ms2, metadata
    def process_peaklist(self, ms1, ms2, metadata):

        self._load_peak_list()
        ms1 = sorted(ms1,key = lambda x: x.mz)
        new_ms1_list = []
        new_ms2_list = []
        new_metadata = {}
        # ms1_mz = [x.mz for z in ms1]
        n_peaks_checked = 0

        ## generate a dict (featid_ms1_dict)to store featid: ms1 pair
        ## O(N) complexisity
        ## build a dict (doc_ms1)for doc_name: ms1 pair first
        doc_ms1, featid_ms1_dict = {}, {}
        for el in ms1:
            doc_name = el.name
            doc_ms1[doc_name] = el
        for k,v in metadata.items():
            if self.id_field and (self.id_field.lower() in v):
                featid = v[self.id_field.lower()]
                featid_ms1_dict[featid] = doc_ms1[k]
            # else:
            #     print(self.id_field)
            #     print(v)

        ## build ms1_ms2 dict, to make searching O(1) in the following loop
        ## key: ms1 object
        ## value: list of ms2
        ms1_ms2_dict = {}
        for el in ms2:
            ms1_ms2_dict.setdefault(el[3], [])
            ms1_ms2_dict[el[3]].append(el)

        if self.id_field and self.csv_id_col: # if the IDs are provided, we match by that
            print("IDs provided ({},{}), using them to match".format(self.id_field,self.csv_id_col))
            match_by_id = True
        else:
            print("IDs not provided, matching on m/z, rt")
            match_by_id = False

        print("Matching peaks...")
        for n_peaks_checked,peak in enumerate(self.ms1_peaks):
            
            if n_peaks_checked % 500 == 0:
                print(n_peaks_checked)
            featid = peak[0]
            peak_mz = peak[1]
            peak_rt = peak[2]
            peak_intensity = None if self.separator in peak[3] else float(peak[3])
            user_cols = peak[4]

            ## first check FeatureId matching
            ## if featureId not exist, then do "mz/rt matching"
            old_ms1 = None

            if match_by_id:
                if featid != None and featid in featid_ms1_dict:
                    old_ms1 = featid_ms1_dict[featid]
            else:
                min_mz = peak_mz - self.mz_tol*peak_mz/1e6
                max_mz = peak_mz + self.mz_tol*peak_mz/1e6
                min_rt = peak_rt - self.rt_tol
                max_rt = peak_rt + self.rt_tol


                ms1_hits = list(filter(lambda x: x.mz >= min_mz and x.mz <= max_mz and x.rt >= min_rt and x.rt <= max_rt,ms1))


                if len(ms1_hits) == 1:
                    # Found one hit, easy
                    old_ms1 = ms1_hits[0]
                elif len(ms1_hits) > 1:
                    # Find the one with the most intense MS2 peak
                    best_ms1 = None
                    best_intensity = 0.0
                    for frag_peak in ms2:
                        if frag_peak[3] in ms1_hits:
                            if frag_peak[2] > best_intensity:
                                best_intensity = frag_peak[2]
                                best_ms1 = frag_peak[3]
                    old_ms1 = best_ms1

            ## Bug fix:
            ## add these two lines to avoid the case that min_ms2_intensity has been set too high,
            ## then most fragments will be removed, and we cannot find a hit for ms1, which will lead to bug:
            ## AttributeError: 'NoneType' object has no attribute 'id'
            if not old_ms1:
                continue

            from time import time
            # make a new ms1 object
            new_ms1 = MS1(old_ms1.id,peak_mz,peak_rt,peak_intensity,old_ms1.file_name,old_ms1.scan_number)
            new_ms1.name = old_ms1.name
            new_ms1_list.append(new_ms1)
            new_metadata[new_ms1.name] = metadata[old_ms1.name]

            ## record user index columns before "mass" column in peaklist file into metadata
            new_metadata[new_ms1.name]['user_cols'] = zip(self.user_cols_names, user_cols)

            if self.separator in peak[3]:
                # print "process sample", str(peak[0]), str(peak[1])
                tokens = []
                for token in peak[3].split(self.separator):
                    try:
                        token = float(token)
                    except:
                        token = None
                    if token <= 0:
                        token = None
                    tokens.append(token)
                # tokens = [float(token) for token in peak[2].split(self.separator)]
                new_metadata[new_ms1.name]['intensities'] = dict(zip(self.sample_names, tokens))

            # Delete the old one so it can't be picked again - removed this, maybe it's not a good idea?
            # pos = ms1.index(old_ms1)
            # del ms1[pos]

            # Change the reference in the ms2 objects to the new ms1 object

            ## Use a dictionary outside the loop to replace the following method, O(N^2) => O(N)
            # ms2_objects = filter(lambda x: x[3] == old_ms1,ms2)
            ms2_objects = []
            if old_ms1 in ms1_ms2_dict:
                ms2_objects = ms1_ms2_dict[old_ms1]

            for frag_peak in ms2_objects:
                new_frag_peak = (frag_peak[0],peak_rt,frag_peak[2],new_ms1,frag_peak[4],frag_peak[5])
                new_ms2_list.append(new_frag_peak)

        # replace the ms1,ms2 and metadata with the new versions
        print("Before {} documents".format(len(ms1)))
        ms1 = new_ms1_list
        ms2 = new_ms2_list
        metadata = new_metadata
        print("Peaklist filtering results in {} documents".format(len(ms1)))
        return ms1, ms2, metadata


    def filter_ms1_intensity(self,ms1,ms2,min_ms1_intensity = 1e6):
        ## Use filter function to simplify code
        print("Filtering MS1 on intensity")
        ## Sometimes ms1 intensity could be None
        ms1 = list(filter(lambda x: False if x.intensity and x.intensity < min_ms1_intensity else True, ms1))
        print("{} MS1 remaining".format(len(ms1)))
        ms2 = list(filter(lambda x: x[3] in set(ms1), ms2))
        print("{} MS2 remaining".format(len(ms2)))
        return ms1, ms2

    def filter_ms2_intensity(self,ms2, min_ms2_intensity = 1e6):
        print("Filtering MS2 on intensity")
        ms2 = list(filter(lambda x: x[2] >= min_ms2_intensity, ms2))
        print("{} MS2 remaining".format(len(ms2)))
        return ms2

    def filter_ms1(self,ms1,ms2,mz_tol = 0.5,rt_tol = 16):
        print("Filtering MS1 to remove duplicates")
        # Filters the loaded ms1s to reduce the number of times that the same molecule has been fragmented


        # Sort the remaining ones by intensity
        ms1_by_intensity = sorted(ms1,key = lambda x: x.intensity,reverse=True)


        final_ms1_list = []
        final_ms2_list = []
        while True:
            if len(ms1_by_intensity) == 0:
                break
            # Take the highest intensity one, find things within the window and remove them
            current_ms1 = ms1_by_intensity[0]
            final_ms1_list.append(current_ms1)
            del ms1_by_intensity[0]

            current_mz = current_ms1.mz
            mz_err = mz_tol*1.0*current_mz/(1.0*1e6)
            min_mz = current_mz - mz_err
            max_mz = current_mz + mz_err

            min_rt = current_ms1.rt - rt_tol
            max_rt = current_ms1.rt + rt_tol

            # find things inside this region
            hits = list(filter(lambda x: x.mz > min_mz and x.mz < max_mz and x.rt > min_rt and x.rt < max_rt,ms1_by_intensity))
            for hit in hits:
                pos = ms1_by_intensity.index(hit)
                del ms1_by_intensity[pos]


        print("{} MS1 remaining".format(len(final_ms1_list)))
        for m in ms2:
            if m[3] in final_ms1_list:
                final_ms2_list.append(m)

        print("{} MS2 remaining".format(len(final_ms2_list)))
        return final_ms1_list,final_ms2_list

    def process_metadata(self, ms1, metadata):
        filtered_metadata = {}
        for m in ms1:
            if m.name in metadata:
                filtered_metadata[m.name] = metadata[m.name]
        metadata = filtered_metadata

        return metadata


# A class to load mzml files
# Will ultimately be able to do method 3

# This method finds each ms2 spectrum in the file and then looks back at the last ms1 scan to find the most
# intense ms1 peak within plus and minus the isolation window. If nothing is found, no document is created
# If it is found, a document is created

# If a peak list is provided it then tries to match the peaks in the peaklist to the ms1 objects, just
# keeping the ms1 objects that can be matched. The matching is done with plus and minus the mz_tol (ppm)
# and plus and minus the rt_tol

## Refactored Sep 21, 2017
## move __init__, and peaklist processing part to parent class *Loader*
class LoadMZML(Loader):
    def __str__(self):
        return "mzML loader"
    def load_spectra(self,input_set):
        import pymzml
        import bisect

        ms1 = []
        ms2 = []
        metadata = {}
        
        ms2_id = 0
        ms1_id = 0


        for input_file in input_set:
            current_ms1_scan_mz = None
            current_ms1_scan_intensity = None
            current_ms1_scan_rt = None
            run = pymzml.run.Reader(input_file, MS1_Precision=5e-6,
                                    extraAccessions=[('MS:1000016', ['value', 'unitName'])],
                                    obo_version='4.0.1')
            file_name = input_file.split('/')[-1]
            previous_precursor_mz = -10
            previous_ms1 = None

            for nc,spectrum in enumerate(run):
                if spectrum['ms level'] == 1:
                    current_ms1_scan_number = nc
                    current_ms1_scan_rt,units = spectrum.scan_time
                    if units == 'minute':
                        current_ms1_scan_rt *= 60.0
                    if current_ms1_scan_rt < self.min_ms1_rt or current_ms1_scan_rt > self.max_ms1_rt:
                        current_ms1_scan_mz = None
                        current_ms1_scan_intensity = None
                    # Note can sometimes get empty scans at the start. If this happens we should ignore.
                    elif len(spectrum.peaks('raw')) > 0:
                        current_ms1_scan_mz,current_ms1_scan_intensity = zip(*spectrum.peaks('raw'))
                    else:
                        current_ms1_scan_mz = None
                        current_ms1_scan_intensity = None

                    previous_precursor_mz = -10
                    previous_ms1 = None
                elif spectrum['ms level'] == 2:
                    # Check that we have an MS1 scan to refer to. If not, skip this one
                    # this can happen if we have blank MS1 scans. We should never get an MS2 scan after a blank MS1
                    # but better to be safe than sorry!
                    if not current_ms1_scan_mz:
                        continue
                    else:
                        precursor_mz = spectrum.selected_precursors[0]['mz']
                        if abs(precursor_mz-previous_precursor_mz) < self.repeated_precursor_match:
                            # Another collision energy perhaps??
                            # if this is the case, we don't bother looking for a parent, but add to the previous one
                            # Make the ms2 objects:
                            if previous_ms1:
                                for mz,intensity in spectrum.centroidedPeaks:
                                    ms2.append((mz,current_ms1_scan_rt,intensity,previous_ms1,file_name,float(ms2_id),nc))
                                    ms2_id += 1
                            else:
                                pass
                        else:
                            # This is a new fragmentation

                            # This finds the insertion position for the precursor mz (i.e. the position one to the right
                            # of the first element it is greater than)
                            precursor_index_ish = bisect.bisect_right(current_ms1_scan_mz,precursor_mz)
                            pos = precursor_index_ish - 1 # pos is now the largest value smaller than ours

                            # Move left and right within the precursor window and pick the most intense parent_scan_number
                            max_intensity = 0.0
                            max_intensity_pos = None
                            while abs(precursor_mz - current_ms1_scan_mz[pos]) < self.isolation_window:
                                if current_ms1_scan_intensity[pos] >= max_intensity:
                                    max_intensity = current_ms1_scan_intensity[pos]
                                    max_intensity_pos = pos
                                pos -= 1
                                if pos < 0:
                                    break
                            pos = precursor_index_ish
                            if pos < len(current_ms1_scan_mz):
                                while abs(precursor_mz - current_ms1_scan_mz[pos]) < self.isolation_window:
                                    if current_ms1_scan_intensity[pos] >= max_intensity:
                                        max_intensity = current_ms1_scan_intensity[pos]
                                        max_intensity_pos = pos
                                    pos += 1
                                    if pos > len(current_ms1_scan_mz)-1:
                                        break
                                # print current_ms1_scan_mz[max_intensity_pos],current_ms1_scan_rt
                            # Make the new MS1 object
                            if (max_intensity > self.min_ms1_intensity) and (not max_intensity_pos == None):
                            # mz,rt,intensity,file_name,scan_number = None):
                                # fix the charge for better loss computation
                                str_charge = spectrum.selected_precursors[0].get('charge',"+1")
                                int_charge = self._interpret_charge(str_charge)

                                # precursormass = current_ms1_scan_mz[max_intensity_pos]
                                parent_mass,single_charge_precursor_mass = self._ion_masses(precursor_mz,int_charge)


                                new_ms1 = MS1(ms1_id,precursor_mz,
                                              current_ms1_scan_rt,max_intensity,file_name,
                                              scan_number = current_ms1_scan_number,
                                              single_charge_precursor_mass = single_charge_precursor_mass)



                                # ms1.append(new_ms1)
                                # ms1_id += 1

                                # Make the ms2 objects:
                                n_found = 0
                                for mz,intensity in spectrum.centroidedPeaks:
                                    if intensity > self.min_ms2_intensity:
                                        ms2.append((mz,current_ms1_scan_rt,intensity,new_ms1,file_name,float(ms2_id),nc))
                                        ms2_id += 1
                                        n_found += 1
                                if n_found > 0:
                                    ms1.append(new_ms1)
                                    ms1_id += 1
                                    metadata[new_ms1.name] = {'most_intense_precursor_mass':current_ms1_scan_mz[max_intensity_pos],
                                                              'parentrt':current_ms1_scan_rt,'scan_number':current_ms1_scan_number,
                                                              'precursor_mass':precursor_mz,'file':file_name,'charge':int_charge,
                                                              'parentmass': parent_mass}


                                    previous_ms1 = new_ms1 # used for merging energies
                                    previous_precursor_mz = new_ms1.mz
               


        print("Found {} ms2 spectra, and {} individual ms2 objects".format(len(ms1),len(ms2)))

        # if self.min_ms1_intensity>0.0:
        #     ms1,ms2 = filter_ms1_intensity(ms1,ms2,min_ms1_intensity = self.min_ms1_intensity)

        # if self.min_ms2_intensity > 0.0:
        #     ms2 = filter_ms2_intensity(ms2, min_ms2_intensity = self.min_ms2_intensity)
        #     # make sure that we haven't ended up with ms1 objects without any ms2
        #     ms1 = []
        #     for m in ms2:
        #         ms1.append(m[3])
        #     ms1 = list(set(ms1))


        if self.peaklist:
            ms1, ms2, metadata = self.process_peaklist(ms1, ms2, metadata)

        if self.duplicate_filter:
            ms1,ms2 = self.filter_ms1(ms1,ms2,mz_tol = self.duplicate_filter_mz_tol,rt_tol = self.duplicate_filter_rt_tol)



        ## class refactor, put filtering inside of the class
        # ms1 = filter(lambda x: x.rt > self.min_ms1_rt and x.rt < self.max_ms1_rt, ms1)
        # ms2 = filter(lambda x: x[3] in set(ms1),ms2)
        # ms2 = filter(lambda x: x[3].rt > self.min_ms1_rt and x[3].rt < self.max_ms1_rt, ms2)

        # Chop out filtered docs from metadata
        metadata = self.process_metadata(ms1, metadata)

        return ms1,ms2,metadata







