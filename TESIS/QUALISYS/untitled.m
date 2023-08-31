clc 
clear
file = readtable("RW001_3D.tsv","FileType","delimitedtext");

file.Properties.VariableNames(1) = "Upstream1_x";
file.Properties.VariableNames(2) = "Upstream1_y";
file.Properties.VariableNames(3) = "Upstream1_z";

file.Properties.VariableNames(4) = "Upstream2_x";
file.Properties.VariableNames(5) = "Upstream2_y";
file.Properties.VariableNames(6) = "Upstream2_z";

file.Properties.VariableNames(7) = "Upstream3_x";
file.Properties.VariableNames(8) = "Upstream3_y";
file.Properties.VariableNames(9) = "Upstream3_z";

file.Properties.VariableNames(10) = "Upstream4_x";
file.Properties.VariableNames(11) = "Upstream4_y";
file.Properties.VariableNames(12) = "Upstream4_z";

file.Properties.VariableNames(13) = "Downstream1_x";
file.Properties.VariableNames(14) = "Downstream1_y";
file.Properties.VariableNames(15) = "Downstream1_z";

file.Properties.VariableNames(16) = "Downstream2_x";
file.Properties.VariableNames(17) = "Downstream2_y";
file.Properties.VariableNames(18) = "Downstream2_z";

file.Properties.VariableNames(19) = "Downstream3_x";
file.Properties.VariableNames(20) = "Downstream3_y";
file.Properties.VariableNames(21) = "Downstream3_z";

file.Properties.VariableNames(22) = "Downstream4_x";
file.Properties.VariableNames(23) = "Downstream4_y";
file.Properties.VariableNames(24) = "Downstream4_z";