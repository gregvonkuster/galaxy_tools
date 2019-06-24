<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            Colony of sample with affy id ${affy_id}
        </h3>
        <table align="center" class="colored">
            %if len(colonies) == 0:
                <tr>
                    <td colspan="2">
                        This sample has no colony
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Latitude</td>
                    <td>Longitude</td>
                    <td>Depth Name</td>
                    <td>Reef ID</td>
                </tr>
                <% ctr = 0 %>
                %for colony in colonies:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${colony[0]}</td>
                        <td>${colony[1]}</td>
                        <td>${colony[2]}</td>
                        <td><a href="${h.url_for( controller='reefs', action='of_colony', sort_id='default', order='default', latitude=colony[0], longitude=colony[1], depth=colony[2], reef_id=colony[3] )}">${colony[3]}</a></td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>
