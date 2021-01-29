function [a, gr, hub] =faceinfo(face)
% Input 'face' - the face id 1 ~ 16
% Output
% 'a' the ranks in competence and popularity dimensions
% 'gr'1,for group1 2 for group2
% 'hub'whether this face plays a role of hub?

        gr2=[2,3,5,8,10,11,13,16]; %gr1=[1,4,6,7,9,12,14,15];
        idxmx=[1 1 1; 2 1 2; 3 1 3; 4 1 4; 5 2 1; 6 2 2; 7 2 3; 8 2 4; 9 3 1; 10 3 2; 11 3 3; 12 3 4; 13 4 1; 14 4 2; 15 4 3; 16 4 4];
        hbc=[7,8,9,10];
        hbp=[2,6,11,15];
        
if numel(face)==1
        gr=ismember(face,gr2)+1;
        a(1) = idxmx(idxmx(:,1)==face,2);
        a(2) = idxmx(idxmx(:,1)==face,3);
        hub(1)=ismember(face,hbc);
        hub(2)=ismember(face,hbp);
    else
        for i=1:numel(face)
            gr(i)=ismember(face(i),gr2)+1;
            a(i,1) = idxmx(idxmx(:,1)==face(i),2);
            a(i,2) = idxmx(idxmx(:,1)==face(i),3);
            hub(i,1)=ismember(face(i),hbc);
            hub(i,2)=ismember(face(i),hbp);
        end
    end 
end
