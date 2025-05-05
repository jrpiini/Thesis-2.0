function LOS = LOS_from_RADec(RA, Dec)

if isscalar(RA)

LOS = [cosd(Dec)*cosd(RA); cosd(Dec)*sind(RA); sind(Dec)];
LOS = LOS / norm(LOS);

else
    LOS = zeros(height(RA), 3);
    for i = 1:height(RA)
        LOS(i,:) = [cosd(Dec(i))*cosd(RA(i)) cosd(Dec(i))*sind(RA(i)) sind(Dec(i))];
        LOS(i,:) = LOS(i,:) / norm(LOS(i,:));
    end
end

end

