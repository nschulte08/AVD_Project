% function to save figures to file:
function fig_save(folder_name, figure_name)
current_folder = pwd;
if ~exist(folder_name,'dir')
    mkdir(folder_name)
end
Figure_folder = sprintf('%s/%s/', current_folder, folder_name);

fig_name =  sprintf('%s.fig',figure_name);
jpg_name =  sprintf('%s.jpg',figure_name);

cd(Figure_folder)
    savefig(fig_name);
    saveas(gcf,jpg_name)
    close(gcf)
cd(current_folder);
end
