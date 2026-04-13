function [system2] = Extract_GSD(gsdname,varargin)

%%NEED to add dihedrals, impropers, and special pairs
%Reads a HOOMD GSD file and writes to a variable "system" with fields
%representing the system. Topology is read from the first timestep in the
%GSD and will not change. 
%
%Optional flags inlcude 'VelocityFlag','DiameterFlag',
%'moment_inertiaFlag', 'ImageFlag', 'OrientationFlag', 'BodyFlag',
%'MassFlag', 'StepFlag'. Deafult will read Stepflag and mass flag, but no other optional
%flags. Set 'AllFlag' to read all components. Set 'VisFlag' to read Mass
%and Diameter. 
%
%Steps to be read can be optionally set as a vector using 'Steps'. Default
%will read all stesp (set to 0). If too many steps are read it will stop at
%the last step. Negative numbers will indicate distance from last step. -1
%is the last step, -2 is the step before, etc. Will only extract 1 stp
%using negative numbers
%The VIADCD options is faster for systems with many timesteps and will
%write a temporary xml and dcd file (using vmd for the dcd file) and then
%load from there.

% py.sys.setdlopenflags(int32(10));
%py.os.chdir(pwd)


p=inputParser;

defaulUnwrappedFlag=0;
defaultAllFlag=0;
defaultVelocityFlag=0;
defaultDiameterFlag=0;
defaultmoment_inertiaFlag=0;
defaultImageFlag=0;
defaultOrientationFlag=0;
defaultBodyFlag=0;
defaultMassFlag=1;
defaultVisFlag=0;
defaultStepFlag=1;
defaultSteps=0;
defaultChargeFlag=0;
defaultDimensionsDynamic=0;
defaultVIADCD=0;
defaultVIADCDframeskip=1;
defaultVIADCDsavelocation='.';
defaultVIADCDfirststep=0;
defaultVIADCDlaststep=-1;

% addParameter(p,'UnwrappedFlag',defaulUnwrappedFlag,@isnumeric)
addParameter(p,'AllFlag',defaultAllFlag,@isnumeric)
addParameter(p,'VelocityFlag',defaultVelocityFlag,@isnumeric)
addParameter(p,'DiameterFlag',defaultDiameterFlag,@isnumeric)
addParameter(p,'moment_inertiaFlag',defaultmoment_inertiaFlag,@isnumeric)
addParameter(p,'ImageFlag',defaultImageFlag,@isnumeric)
addParameter(p,'OrientationFlag',defaultOrientationFlag,@isnumeric)
addParameter(p,'BodyFlag',defaultBodyFlag,@isnumeric)
addParameter(p,'MassFlag',defaultMassFlag,@isnumeric)
addParameter(p,'VisFlag',defaultVisFlag,@isnumeric)
addParameter(p,'StepFlag',defaultStepFlag,@isnumeric)
addParameter(p,'Steps',defaultSteps,@isnumeric)   
addParameter(p,'ChargeFlag',defaultChargeFlag,@isnumeric)
addParameter(p,'DimensionsDynamic',defaultDimensionsDynamic,@isnumeric)
addParameter(p,'VIADCD',defaultVIADCD,@isnumeric)
addParameter(p,'VIADCDframeskip',defaultVIADCDframeskip,@isnumeric)
addParameter(p,'VIADCDsavelocation',defaultVIADCDsavelocation,@ischar)
addParameter(p,'VIADCDfirststep',defaultVIADCDfirststep,@isnumeric)
addParameter(p,'VIADCDlaststep',defaultVIADCDlaststep,@isnumeric)
parse(p,varargin{:});

fields=fieldnames(p.Results);
for i=1:numel(fields)
    eval([fields{i} '=' 'p.Results.' fields{i} ';']);
end
    


if AllFlag==1
%     UnwrappedFlag=1;
    VelocityFlag=1;
    DiameterFlag=1;
    moment_inertiaFlag=1;
    ImageFlag=1;
    OrientationFlag=1;
    BodyFlag=1;
    MassFlag=1;
    StepFlag=1;
    ChargeFlag=1;
    if VisFlag==1
        error("Vis and All Flag can't both be set");
    end
end

if VisFlag==1
%     UnwrappedFlag=1;
%     VelocityFlag=1;
    DiameterFlag=1;
%     MOIFlag=1;
%     ImageFlag=1;
%     OrientationFlag=1;
%     BodyFlag=1;
    MassFlag=1;
    if AllFlag==1
        error("Vis and All Flag can't both be set");
    end
end


if VIADCD==1
    %writesavedcd.tcl
    if ~endsWith(VIADCDsavelocation,'/')
        VIADCDsavelocation=[VIADCDsavelocation '/']
    end
    fidVIADCD=fopen('savedcd.tcl','wt')
    fprintf(fidVIADCD,['mol new {' gsdname '} type {gsd} first ' num2str(VIADCDfirststep) ' last ' num2str(VIADCDlaststep) ' step ' num2str(VIADCDframeskip) ' waitfor -1\nset id [molinfo top get id ]\nset numframes [molinfo $id get numframes]\nset numframes [expr $numframes - 1]\nanimate write dcd {' VIADCDsavelocation 'VIADCD.dcd} beg 0 end $numframes skip 0 0'])
    fclose(fidVIADCD)
    %system(['vmd -dispdev text ' gsdname ' -eofexit < ./GIT/matlab_gsd_reader/savedcd.tcl > output.log'])
%     system(['vmd -dispdev text -eofexit < savedcd.tcl > output.log'])
    system(['vmd -dispdev text -eofexit < savedcd.tcl'])
    system3=Extract_GSD(gsdname,'Steps',1);
    Write_XML(0,[VIADCDsavelocation 'VIADCD.xml'],system3);
    system2=Extract_XML_DCD([VIADCDsavelocation 'VIADCD.xml'],[VIADCDsavelocation 'VIADCD.dcd'],0,Steps);
    
    delete([VIADCDsavelocation 'VIADCD.xml']);
    delete([VIADCDsavelocation 'VIADCD.dcd']);
%     delete('savedcd.tcl');
elseif VIADCD ~=0
    error('VIADCD wrongly specified')
else







% open GSD file
t = py.gsd.hoomd.open(gsdname);

% step 1 — get a handle to Python's __getitem__ method
getitem = py.getattr(t, '__getitem__');

% step 2 — call that method with argument 0
s = getitem(int32(0));

% now s is a gsd.hoomd.Frame
system2.natoms = double(s.particles.N);



   % MATLAB → Python indexing, frame 0

    %s=t.read_frame(py.len(t)-1);

    %natoms
    system2.natoms=s.particles.N+1-1;

    %attype
    for i=1:length(s.particles.types)
        a=char(s.particles.types(i));
        types(i,:)=a(3:end-2);
    end
    % types=cellfun(@char,cell(s.particles.types));
    attype = double(py.array.array('d',py.numpy.nditer(s.particles.typeid)));
    system2.attype=types(attype+1,:);

    %bonds
    if s.bonds.N+1-1~=0
        system2.bonds = double(py.array.array('d',py.numpy.nditer(s.bonds.group)));
        system2.bonds =reshape(system2.bonds.' ,[2 s.bonds.N+1-1])';
    else
        system2.bonds = [];
    end

    %dimensions
    system2.dim = double(py.array.array('d',py.numpy.nditer(s.configuration.box)));
    system2.dim=system2.dim(1:3);

    %mass
    if MassFlag==1
        system2.mass = double(py.array.array('d',py.numpy.nditer(s.particles.mass)))';
    end

    %angles
    if s.angles.N+1-1~=0
        system2.angles = double(py.array.array('d',py.numpy.nditer(s.angles.group)));
        system2.angles=reshape(system2.angles.' ,[3 s.angles.N+1-1])';
    else
        system2.angles=[];
    end

    %Diameter 
    if DiameterFlag==1
        system2.diameter = double(py.array.array('d',py.numpy.nditer(s.particles.diameter)))';
    end

    %Charge
    if ChargeFlag==1
        system2.charge = double(py.array.array('d',py.numpy.nditer(s.particles.charge)))';
    end

    %body
    if BodyFlag==1
        system2.body = double(py.array.array('d',py.numpy.nditer(s.particles.body)))';
    end


     %dynamic data
     if Steps==0
         stepslist=0:py.len(t)-1;
     elseif max(Steps)>py.len(t)-1
         stepslist=Steps(find(Steps<=py.len(t)-1));
     elseif Steps<0
         stepslist=py.len(t)+Steps; %if negative go that many steps from the end
     else
         stepslist=Steps;
     end
     stepindex=0;

     %stepsize
getitem = py.getattr(t, '__getitem__');  % Python [] operator

if py.len(t) > 1
    % --- Read first two frames to compute step size ---
    s1 = getitem(int32(0));
    ts1 = double(s1.configuration.step);

    s2 = getitem(int32(1));
    ts2 = double(s2.configuration.step);

    if length(stepslist) > 1
        system2.stepsize = (ts2 - ts1) * (stepslist(2) - stepslist(1));
    else
        system2.stepsize = NaN;
    end

    system2.initialstep = ts1;
    
    % --- Read last frame to get current timestep ---
    s_last = getitem(int32(stepslist(end)));
    system2.timestep = double(s_last.configuration.step);

else
    % --- Only one frame in trajectory ---
    s1 = getitem(int32(0));
    system2.timestep = double(s1.configuration.step);
    system2.stepsize = NaN;
end


    for i = stepslist
    stepnum = i;
    if mod(stepnum, 1000) == 0 && stepnum ~= 0
        disp(stepnum);
    end
    stepindex = stepindex + 1;

    % --- Get frame i ---
    s = getitem(int32(i));

    % --- Extract positions ---
    % Convert numpy array → MATLAB double array
    pos_np = s.particles.position;
    temppos = double(py.array.array('d', py.numpy.nditer(pos_np)))';
    system2.pos(:, :, stepindex) = reshape(temppos, [3, system2.natoms])';

    % --- Extract velocities (if present) ---
    if VelocityFlag == 1
        vel_np = s.particles.velocity;
        tempvel = double(py.array.array('d', py.numpy.nditer(vel_np)))';
        system2.vel(:, :, stepindex) = reshape(tempvel, [3, system2.natoms])';
    end
end

        if ImageFlag==1
            tempimg = double(py.array.array('d',py.numpy.nditer(s.particles.image)))';
            system2.img(:,:,stepindex)= reshape(tempimg,[3,system2.natoms])';
        end

        if moment_inertiaFlag==1
            tempMOI = double(py.array.array('d',py.numpy.nditer(s.particles.moment_inertia)))';
            system2.moment_inertia(:,:,stepindex)= reshape(tempMOI,[3,system2.natoms])';
        end

        if OrientationFlag==1
            tempO = double(py.array.array('d',py.numpy.nditer(s.particles.orientation)))';
            system2.orientation(:,:,stepindex)= reshape(tempO,[4,system2.natoms])';
        end
        if DimensionsDynamic==1
            dimtemp=double(py.array.array('d',py.numpy.nditer(s.configuration.box)));
            system2.dim(:,:,stepindex)=dimtemp(1:3);
        end
    end
 

end


 

