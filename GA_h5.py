import h5py
import numpy as np
import os

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

                        h5group.create_dataset(key, data=val)
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
    constrains = f"g{global_params['switch_to_exc_drive_for_dsi']}_i{global_params['RF_constrains']['excitation']['inactivation']}_s{global_params['RF_constrains']['excitation']['doSurround']}_keepTau_{global_params['RF_constrains']['excitation']['same cs tau']}_"
    constrains += f'c_amp{global_params["RF_constrains"]["excitation"]["center"]["varyAmplitude"]}_c_kin{global_params["RF_constrains"]["excitation"]["center"]["varyKinetics"]}_c_size{global_params["RF_constrains"]["excitation"]["center"]["varySize"]}_c_or{global_params["RF_constrains"]["excitation"]["center"]["varyOrientation"]}'
    constrains += f's_amp{global_params["RF_constrains"]["excitation"]["surround"]["varyAmplitude"]}_s_kin{global_params["RF_constrains"]["excitation"]["surround"]["varyKinetics"]}_s_size{global_params["RF_constrains"]["excitation"]["surround"]["varySize"]}_s_or{global_params["RF_constrains"]["excitation"]["surround"]["varyOrientation"]}'
    
    if(global_params['cellType'] in ['SAC','SAC network']):
        file_name= f"_s{global_params['numSpeed']}_c{global_params['numContrast']}_d{global_params['numDir']}_DIST{global_params['distSAC']}_NET{global_params['numSAClayersX']}x{global_params['numSAClayersY']}_{global_params['activeChannels']}_{global_params['network']}_j{global_params['job_id']}.h5"
    else:
        file_name= f"_s{global_params['numSpeed']}_{constrains}_j{global_params['job_id']}.h5"
    #if global_params['use_exc_drive_for_dsi']:

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
        models[0].output_params['cell']['somaV_all_angles']= pad_array(global_params, 0, models[0].output_params['cell']['somaV_all_angles'])

        # save presynaptic inputs
        for bc in range(len(models[0].output_params['input_exc_cell'])):
            models[0].output_params['input_exc_cell'][bc]['somaV_all_angles']= pad_array(global_params, -1, models[0].output_params['input_exc_cell'][bc]['somaV_all_angles'])
            models[0].output_params['input_exc_cell'][bc]['cell_input_all_angles']= pad_array(global_params, -1, models[0].output_params['input_exc_cell'][bc]['cell_input_all_angles'])
            models[0].output_params['input_exc_cell'][bc]['cell_output_all_angles']= pad_array(global_params, -1, models[0].output_params['input_exc_cell'][bc]['cell_output_all_angles'])

        for ac in range(len(models[0].output_params['input_inh_cell'])):
            models[0].output_params['input_inh_cell'][ac]['somaV_all_angles']= pad_array(global_params, -1, models[0].output_params['input_inh_cell'][ac]['somaV_all_angles'])
            models[0].output_params['input_inh_cell'][ac]['cell_input_all_angles']= pad_array(global_params, -1, models[0].output_params['input_inh_cell'][ac]['cell_input_all_angles'])
            models[0].output_params['input_inh_cell'][ac]['cell_output_all_angles']= pad_array(global_params, -1, models[0].output_params['input_inh_cell'][ac]['cell_output_all_angles'])


        # Remove large structures
        if not global_params['save_all']:
            models[0].output_params['trajectory']= []
            models[0].input_params['trajectory']= []
            models[0].input_params['zero_trajectory']= []
            
             
        output_grp = f.create_group(f"output")
        save_to_group(output_grp, models[0].output_params)
            
        
        # Inputs   
        input_grp = f.create_group(f"input")
        save_to_group(input_grp, models[0].input_params)
        
        # Params
        param_grp = f.create_group(f"params")
        save_to_group(param_grp, global_params)  

        f.close()      

