<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            Reef of Colony
        </h3>
        <h4>
            Latitude: <b>${latitude}</b><br/>
            Longitude: <b>${longitude}</b><br/>
            Depth: <b>${depth}</b><br/>
        </h4>
        <table align="center" class="colored">
            %if len(reefs) == 0:
                <tr>
                    <td colspan="2">
                        There is no reef for this colony
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Name</td>
                    <td>Region</td>
                    <td>Latitude</td>
                    <td>Longitude</td>
                    <td>GPS Coordinates Associated With</td>
                </tr>
                <% ctr = 0 %>
                %for reef in reefs:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${reef[0]}</td>
                        <td>${reef[1]}</td>
                        <td>${reef[2]}</td>
                        <td>${reef[3]}</td>
                        <td>${reef[4]}</td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>
