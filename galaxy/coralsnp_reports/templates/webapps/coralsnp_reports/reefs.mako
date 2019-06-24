<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            All reefs for uploaded samples
        </h3>
        <table align="center" class="colored">
            %if len(reefs) == 0:
                <tr>
                    <td colspan="2">
                        There are no reefs for uploaded samples
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Name</td>
                    <td>Region</td>
                    <td>Latitude</td>
                    <td>Longitude</td>
                    <td>Geographic Origin</td>
                    <td>Samples of this Reef</td>
                </tr>
                <% ctr = 0 %>
                %for reef in reefs:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${reef[1]}</td>
                        <td>${reef[2]}</td>
                        <td>${reef[3]}</td>
                        <td>${reef[4]}</td>
                        <td>${reef[5]}</td>
                        <td><a href="${h.url_for( controller='samples', action='of_reef', sort_id='default', order='default', reef_id=reef[0], name=reef[1], region=reef[2], latitude=reef[3], longitude=reef[4], geographic_origin=reef[5] )}">Samples</a></td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

