function computer_name=machine
% computer_name = machine;
%
% This function returns the machine name for a PC running
% Win95/98, NT4.0 or Win 2000. Output is a string of 
% form '\\mycomputer'.
%
% (c) Dr. Ben Quine, Dept. of Physics,
% University of Toronto, 23-MAY-2000.
%
% Email: ben@atmosp.physics.utoronto.ca
%

if isempty(findstr(computer,'PCWIN'))
   error('Not a PC Windows machine. Cannot Determine computer name.')
end   
% finding the computer name (win 95/98)
[status, computer_name]=dos('net config');
nameind = findstr(computer_name,'\\');
if isempty(nameind)
  % finding the computer name (win NT4.0, win 2000)
	[status, computer_name]=dos('net config workstation');
	nameind = findstr(computer_name,'\\');
end	
if isempty(nameind)
   error('Cannot Determine computer name.')
end
computer_name = computer_name(nameind:length(computer_name));
nameind = findstr(computer_name,char(10));
% record computer name
computer_name = computer_name(1:(nameind-1));
