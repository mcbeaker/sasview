
; Script generated by the Inno Setup Script Wizard

; and local_config.py located in this directory.
 ; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!
[Setup]

ChangesAssociations = yes
AppName = SansView-Dev04192011
AppVerName = SansView-Dev04192011
AppPublisher = (c) 2009, University of Tennessee
AppPublisherURL = http://danse.chem.utk.edu
AppSupportURL = http://danse.chem.utk.edu
AppUpdatesURL = http://danse.chem.utk.edu 
DefaultDirName = {pf}\SansView-Dev04192011
DefaultGroupName = DANSE\SansView-Dev04192011
DisableProgramGroupPage = yes
LicenseFile = license.txt
OutputBaseFilename = setupSansView
SetupIconFile = images\ball.ico
Compression = lzma
SolidCompression = yes
PrivilegesRequired = none


[Registry]
Root: HKCR;	Subkey: ".xml";	ValueType: string;	ValueName: "";	ValueData: "{app}\SansView.exe";	 Flags: uninsdeletevalue
Root: HKCR;	Subkey: ".asc";	ValueType: string;	ValueName: "";	ValueData: "{app}\SansView.exe";	 Flags: uninsdeletevalue
Root: HKCR;	Subkey: ".dat";	ValueType: string;	ValueName: "";	ValueData: "{app}\SansView.exe";	 Flags: uninsdeletevalue
Root: HKCR;	Subkey: ".tif";	ValueType: string;	ValueName: "";	ValueData: "{app}\SansView.exe";	 Flags: uninsdeletevalue
Root: HKCR;	Subkey: ".tiff";	ValueType: string;	ValueName: "";	ValueData: "{app}\SansView.exe";	 Flags: uninsdeletevalue
Root: HKCR;	Subkey: ".sans";	ValueType: string;	ValueName: "";	ValueData: "{app}\SansView.exe";	 Flags: uninsdeletevalue
Root: HKCR;	Subkey: ".svs";	ValueType: string;	ValueName: "";	ValueData: "{app}\SansView.exe";	 Flags: uninsdeletevalue
Root: HKCR;	Subkey: ".fitv";	ValueType: string;	ValueName: "";	ValueData: "{app}\SansView.exe";	 Flags: uninsdeletevalue
Root: HKCR;	Subkey: ".inv";	ValueType: string;	ValueName: "";	ValueData: "{app}\SansView.exe";	 Flags: uninsdeletevalue
Root: HKCR;	Subkey: ".prv";	ValueType: string;	ValueName: "";	ValueData: "{app}\SansView.exe";	 Flags: uninsdeletevalue
Root: HKCR; Subkey: "{app}\images\ball.ico";	ValueType: string; ValueName: "";	ValueData: "{app}\SansView.exe,0"
Root: HKCR; Subkey: "{app}\SansView.exe\shell\open\command";	ValueType: string; ValueName: "";	ValueData: """{app}\SansView.exe""  ""%1"""


[Languages]
Name: "english";	MessagesFile: "compiler:Default.isl"


[Tasks]
Name: "desktopicon";	Description: "{cm:CreateDesktopIcon}";	GroupDescription: "{cm:AdditionalIcons}";	Flags: unchecked


[Files]
Source: "dist\SansView.exe";	DestDir: "{app}";	Flags: ignoreversion
Source: "dist\*";	DestDir: "{app}";	Flags: ignoreversion recursesubdirs createallsubdirs
Source: "images\*";	DestDir: "{app}\images";	Flags: ignoreversion recursesubdirs createallsubdirs
Source: "test\*";	DestDir: "{app}\test";	Flags: ignoreversion recursesubdirs createallsubdirs
;	NOTE: Don't use "Flags: ignoreversion" on any shared system files

[Icons]
Name: "{group}\SansView-Dev04192011";	Filename: "{app}\SansView.exe";	WorkingDir: "{app}" 
Name: "{group}\{cm:UninstallProgram, SansView-Dev04192011}";	 Filename: "{uninstallexe}" 
Name: "{commondesktop}\SansView-Dev04192011";	Filename: "{app}\SansView.exe";	Tasks: desktopicon; WorkingDir: "{app}" 


[Run]
Filename: "{app}\SansView.exe";	Description: "{cm:LaunchProgram, SansView-Dev04192011}";	Flags: nowait postinstall skipifsilent
