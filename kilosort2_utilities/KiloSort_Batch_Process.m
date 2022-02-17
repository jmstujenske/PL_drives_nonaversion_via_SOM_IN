%Create a binary file and run kilosort on all the files

file_dir='C:\Users\Admin\Dropbox (ListonLab)\DataForJoe_SOMcreBLArec\';
data_dir='DFC Protocol\';
animal_dirs=dir([file_dir data_dir]);
animal_dirs=animal_dirs(3:end);
day_headers={'HAB','PreHabLaser','PreRecallLaser','RECALL','TRAINING'};
exclude_electrodes=[32];
for animal_rep=1:length(animal_dirs)
    for day_rep=4
        clear mat_file
        day_dir=dir([file_dir data_dir animal_dirs(animal_rep).name,'\','*',day_headers{day_rep}]); 
        day_dir=[day_dir.folder,'\',day_dir.name];
        test_bin=dir([day_dir,'\*.bin']);
        if isempty(test_bin)
            
        mat_file=dir([day_dir,'\*.mat']);
        if ~isempty(mat_file)
            size_files=vertcat(mat_file.bytes);
            if length(size_files)>1
                [~,ind]=max(size_files);
            else
                ind=1;
            end
        mat_file=[mat_file(ind).folder,'\',mat_file(ind).name];
        animal_data=load(mat_file);
        electrodenumbers=vertcat(animal_data.NS5.ElectrodesInfo.ElectrodeID);
        
        fid = fopen([day_dir,'\data_binary.bin'],'w');
        fwrite(fid, int16(animal_data.NS5.Data(~ismember(electrodenumbers,exclude_electrodes),:)) ,'int16');
        fclose(fid);
        end
        end
        test_npy=dir([day_dir,'\*.npy']);
        test_bin=dir([day_dir,'\*.bin']);
        if isempty(test_npy) && ~isempty(test_bin)
%         if ~isempty(test_bin)
            run_kilosort_onfile(day_dir);
        end
end
end