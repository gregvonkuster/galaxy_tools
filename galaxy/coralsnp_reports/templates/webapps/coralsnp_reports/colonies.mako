<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            All colonies for uploaded samples
        </h3>
        <table align="center" class="colored">
            %if len(colonies) == 0:
                <tr>
                    <td colspan="2">
                        There are no colonies for uploaded samples
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Latitude</td>
                    <td>Longitude</td>
                    <td>Depth</td>
                    <td>Reef ID</td>
                    <td>Samples of this Colony</td>
                </tr>
                <% ctr = 0 %>
                %for colony in colonies:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${colony[1]}</td>
                        <td>${colony[2]}</td>
                        <td>${colony[3]}</td>
                        <td><a href="${h.url_for( controller='reefs', action='of_colony', sort_id='default', order='default', latitude=colony[1], longitude=colony[2], depth=colony[3], reef_id=colony[4] )}">${colony[4]}</a></td>
                        <td><a href="${h.url_for( controller='samples', action='of_colony', sort_id='default', order='default', colony_id=colony[0], latitude=colony[1], longitude=colony[2], depth=colony[3], reef_id=colony[4] )}">Samples</a></td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

