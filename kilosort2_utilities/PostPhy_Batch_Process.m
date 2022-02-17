%%After cultivating in Phy, run this script to extract all of the event
%%time stamps and spike data

file_dir='C:\Users\Admin\XXXXXXX\';
data_dir='XXXXXXX\';
animal_dirs=dir([file_dir data_dir]);
animal_dirs=animal_dirs(3:end);
day_headers={'RECALL'};
exclude_electrodes=[];%specify electrodes that were excluded
clear data
for animal_rep=1:n_animal
    for day_rep=4
        clear mat_file
        day_dir=dir([file_dir data_dir animal_dirs(animal_rep).name,'\','*',day_headers{day_rep}]);
        day_name=day_dir.name;
        if exist([seconddir,day_name],'file')
            day_dir=[seconddir,day_name];
        else
            day_dir=[day_dir.folder,'\',day_dir.name];
        end
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
        data(animal_rep)=process_post_phy(day_dir,animal_data,exclude_electrodes);
        end
    end
end
