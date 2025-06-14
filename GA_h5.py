import h5py
import numpy as np
import os

def save_dict_to_hdf5(h5group, data):
    def save_list(h5group, key, val):
            #list_group = h5group.create_group(key)
        if all(isinstance(item, np.ndarray) for item in val):
            print("dict-array",  key)
            for i, item in enumerate(data):
                h5group.create_dataset(key+"_"+str(i), data=item)                    
        else:
            print("dict-not array", key, type(val))
            
            if isinstance(val, list):
                try:
                    if np.array(val).dtype.kind in {'i', 'f'}:
                        print("All numeric (int or float).",key)
                        h5group.create_dataset(key,data=val) 
                   
                    else:
                        print("Contains non-numeric data.",key)
                        save_dict_to_hdf5(h5group, val)
                except Exception:
                    print("Exception",key)
                    save_dict_to_hdf5(h5group, val)    
                
                #for item in val:
                #    print(type(item),isinstance(item, np.ndarray),isinstance(item, int))
                #if all((isinstance(item, np.ndarray) or isinstance(item, int)) for item in val):
                #    print("number")
                #else:
                #    print("not num")

    print("__",str(h5group),type(data))
    
    if isinstance(data, dict):
        print("dict",str(h5group),type(data))
        for key, val in data.items():
            key = str(key)
            if isinstance(val, dict):
                subgroup = h5group.create_group(key)
                save_dict_to_hdf5(subgroup, val)
            elif isinstance(val, list):
                pass
                
                #save_list(h5group, key, val)
                #for i, item in enumerate(val):
                    #item_group = list_group.create_group(str(i))
                 #   if isinstance(item, np.ndarray):
                  #      print(list_group, "array", key)
                   # else:
                    #    print(list_group, "not array", key, item)
                       
                    #save_dict_to_hdf5(list_group, item, key)  
            else:
                print("simple",str(h5group),type(val))
                h5group.create_dataset(key, data=val)              
    elif isinstance(data, list):
        #save_list(h5group, "_", data)
        
        print("LIST",str(h5group),type(data))
        #list_group = h5group.create_group(name)
        #if all(isinstance(item, np.ndarray) for item in data):
        #    print("all_array",str(h5group),type(data))
        for i, item in enumerate(data):
            h5group.create_dataset(str(i), data=item)
        #else:
            #print("not_array_",str(h5group),type(data))
            #for i, item in enumerate(data):                
            #    item_group = list_group.create_group(str(i))
            #    save_dict_to_hdf5(item_group, item, str(i))                  

    #elif isinstance(data, np.ndarray):
    #    print("_array_",str(h5group),str(data))
    #    h5group.create_dataset(name, data=data)
    else:
        #print("simple",str(h5group),str(data))
        h5group.create_dataset("_", data=data)     
"""
                if all(isinstance(item, dict) for item in val):
                    list_group = h5group.create_group(key)
                    for i, item in enumerate(val):
                        item_group = list_group.create_group(str(i))
                        save_dict_to_hdf5(item_group, item)
                elif all(isinstance(item, np.ndarray) for item in val) and len({item.shape for item in val}) > 1:
                    list_group = h5group.create_group(key)
                    for i, item in enumerate(val):
                        list_group.create_dataset(str(i), data=item)
                else:
                    try:
                        arr = np.array(val)
                        h5group.create_dataset(key, data=arr)
                    except Exception:
                        h5group.create_dataset(key, data=str(val))
            elif isinstance(val, (np.ndarray, int, float, str, np.integer, np.floating)):
                h5group.create_dataset(key, data=val)
            else:
                try:
                    arr = np.array(val)
                    h5group.create_dataset(key, data=arr)
                except Exception:
                    h5group.create_dataset(key, data=str(val))
    #elif isinstance(data, list):
    #    print("list", data)
    #    for i, item in enumerate(data):
    #        item_group = h5group.create_group(str(i))
     #       save_dict_to_hdf5(item_group, item)        
    else:
        try:
            if all(isinstance(item, np.ndarray) for item in data) and len({item.shape for item in data}) > 1:
                h5group.create_dataset("value", data=np.array(data))
            else:
                for i, item in enumerate(data):
                    h5group.create_dataset(str(i), data=data)
        except Exception:
                #list_group = h5group.create_group(key)
            
            h5group.create_dataset("value", data=str(data)) 
            """

def load_dict_from_hdf5(h5group):
    result = {}
    for key in h5group:
        item = h5group[key]
        if isinstance(item, h5py.Group):
            if all(k.isdigit() for k in item.keys()):
                # this is a list-like group
                sorted_items = [item[k] for k in sorted(item.keys(), key=int)]
                result[key] = [np.array(ds) for ds in sorted_items]
            else:
                result[key] = load_dict_from_hdf5(item)
        else:
            result[key] = np.array(item)
    return result


def save_class_to_hdf5(filename, obj):
    with h5py.File(filename, 'w') as f:
        save_to_group(f, obj.__dict__)

def save_to_group(h5group, data):
    if isinstance(data, dict):
        for key, val in data.items():
            key = str(key)
            if isinstance(val, dict):
                subgroup = h5group.create_group(key)
                save_to_group(subgroup, val)
            elif isinstance(val, list):
                if all(isinstance(item, dict) for item in val):

                    list_group = h5group.create_group(key)
                    for i, item in enumerate(val):
                        item_group = list_group.create_group(str(i))
                        save_to_group(item_group, item)
                elif all(isinstance(item, np.ndarray) for item in val) or all(isinstance(item, (int, float)) for item in val):
                    try:

                        h5group.create_dataset(key, data=val)#np.array(val, dtype=object), dtype=h5py.special_dtype(vlen=np.float64))
                    except:
                        for i, item in enumerate(val):
                            h5group.create_dataset(f"{key}_{i}", data=item)
                else:

                    h5group.create_dataset(key, data=val)

            elif isinstance(val, (np.ndarray, int, float, str, np.integer, np.floating)):
                h5group.create_dataset(key, data=val)
            else:
                try:
                    h5group.create_dataset(key, data=np.array(val))
                except Exception:
                    h5group.create_dataset(key, data=(val))
    else:
        print("Top-level data must be a dictionary")

def pad_array(global_params, cell, arrays):
    padded=[]
    if((global_params['save_all']) or (cell == 0)):
        if(len(arrays) > 0):
            max_len = max(len(a) for a in arrays)            
            # Pad each array with np.nan to match max_len
            padded = np.full((len(arrays), max_len), np.nan)
            for i, arr in enumerate(arrays):
                padded[i, :len(arr)] = arr       
    return padded
    

def GA_SaveH5(global_params, models, final):
    if(global_params['cellType'] in ['SAC','SAC network']):
        file_name= f"_s{global_params['numSpeed']}_c{global_params['numContrast']}_d{global_params['numDir']}_DIST{global_params['distSAC']}_NET{global_params['numSAClayersX']}x{global_params['numSAClayersY']}_{global_params['activeChannels']}_j{global_params['job_id']}.h5"
    else:
        file_name= f"_s{global_params['numSpeed']}_c{global_params['numContrast']}_d{global_params['numDir']}_{global_params['cellType']}_j{global_params['job_id']}.h5"
    if(final):
        file_start= "result"
        try:
            file_path= f"results/temp{file_name}"
            os.remove(file_path)
            print("Temp file deleted successfully.")
        except FileNotFoundError:
            print("Temp file does not exist.")
        except Exception as e:
            print(f"Error deleting file: {e}")
    else:
        file_start= "temp"
   
    file_full=f"results/{file_start}{file_name}"

    with h5py.File(file_full, "w") as f:


        
        # Outputs, fix uneven lengths
        for cell in range(global_params['numCells']):
            
            models[0].output_params['cell'][cell]['dendrite_vectors']['ca_right_all_angles']= pad_array(global_params, cell, models[0].output_params['cell'][cell]['dendrite_vectors']['ca_right_all_angles'])   
            models[0].output_params['cell'][cell]['dendrite_vectors']['ca_left_all_angles']= pad_array(global_params, cell, models[0].output_params['cell'][cell]['dendrite_vectors']['ca_left_all_angles'])   
            models[0].output_params['cell'][cell]['dendrite_vectors']['v_right_all_angles']= pad_array(global_params, cell, models[0].output_params['cell'][cell]['dendrite_vectors']['v_right_all_angles'])   
            models[0].output_params['cell'][cell]['dendrite_vectors']['v_left_all_angles']= pad_array(global_params, cell, models[0].output_params['cell'][cell]['dendrite_vectors']['v_left_all_angles'])   
            models[0].output_params['cell'][cell]['somaV_all_angles']= pad_array(global_params, cell, models[0].output_params['cell'][cell]['somaV_all_angles'])
            #print(type(models[0].output_params['cell'][cell]['conductance']['k_total']))
            models[0].output_params['cell'][cell]['conductance']['k_total']= pad_array(global_params, cell, models[0].output_params['cell'][cell]['conductance']['k_total'])   
            models[0].output_params['cell'][cell]['conductance']['ca_total']= pad_array(global_params, cell, models[0].output_params['cell'][cell]['conductance']['ca_total'])   
                
        models[0].output_params['InputE_syn_time_copy']= pad_array(global_params, cell, models[0].output_params['InputE_syn_time_copy'])  
        models[0].output_params['InputI_syn_time_copy']= pad_array(global_params, cell, models[0].output_params['InputI_syn_time_copy'])  
        models[0].output_params['I_syn_g']= pad_array(global_params, cell, models[0].output_params['I_syn_g'])  

        # Do not save these
        models[0].output_params['cell'][cell]['conductance']['k_all']= []
        models[0].output_params['cell'][cell]['conductance']['ca_all']= []
              
        output_grp = f.create_group(f"output")
        save_to_group(output_grp, models[0].output_params)
            
        
        # Inputs   
        input_grp = f.create_group(f"input")
        save_to_group(input_grp, models[0].input_params)
        
        # Params
        param_grp = f.create_group(f"params")
        save_to_group(param_grp, global_params)  

        f.close()      

def GA_LoadInputsFromH5(model):
    with h5py.File("inputs.h5", "r") as f:
        model.input_params = load_dict_from_hdf5(f["input"])