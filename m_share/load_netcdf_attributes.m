function [Att,GroupAtt] = load_netcdf_attributes(filename)
% Att = load_netcdf_attributes(filename)
% fetch all global attributes from an NC file to a structure Att.
% Field % names are the attribute names.
% Required input:
%    filename: string of the filename including .nc OR it's the output from
%    ncinfo i.e., a structure.
% Ilker Fer, 20220512
% CBluteau - dealt with attribute names starting with underscores (no idea
% why they are writing them that way!)
%      Added ability to read attributes at group level
%

if isstruct(filename)
    ncInfo=filename;
else
    ncInfo=ncinfo(filename);
end

%%
nGr=length(ncInfo.Groups);

Att=read_atts(ncInfo);
if nGr==0
    GroupAtt=struct([]);
    return;
end

for kk=1:nGr
    A=ncInfo.Groups(kk);
    GroupAtt.(A.Name)=read_atts(A);
end


function Att=read_atts(nInfo)
% nInfo could be the output from e.g., nInfo=ncinfo; or  group
% info=nInfo.Groups(gid) (with gid=1, 2  etc) or even varInfo=nInfo.Variables(vId)
% Outputs all the attributes associated with the passed info structure.
nAtt=length(nInfo.Attributes);
if nAtt==0
    Att=struct([]);
else
    for I=1:nAtt
        name=clean_names_attributes(nInfo.Attributes(I).Name);
        val = nInfo.Attributes(I).Value;

        Att.(name)=val;
        clear name val

    end
end


function newVarName=clean_names_attributes(varName)
% Takes any string, and cleans it up such that it can be used as a
% fieldname in a matlab structure. For instance, :, -, and prefixed _ are
% removed.


tmp = regexprep(varName,'\s','_');
tmp = regexprep(tmp,'^_+','');
newVarName= regexprep(tmp,'-','_');


