function [UTC] = JD_to_UTC(JD)

UTC = [(year(datetime(JD, 'ConvertFrom', 'juliandate'))), (month(datetime(JD, 'ConvertFrom', 'juliandate'))), (day(datetime(JD, 'ConvertFrom', 'juliandate'))), (hour(datetime(JD, 'ConvertFrom', 'juliandate'))), (minute(datetime(JD, 'ConvertFrom', 'juliandate'))), (second(datetime(JD, 'ConvertFrom', 'juliandate')))];

end

