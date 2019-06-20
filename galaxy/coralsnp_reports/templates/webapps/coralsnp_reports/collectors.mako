<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            All sample collectors
        </h3>
        <table align="center" class="colored">
            %if len(collectors) == 0:
                <tr>
                    <td colspan="2">
                        There are no sample collectors
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Last Name</td>
                    <td>First Name</td>
                    <td>Organization</td>
                    <td>Email</td>
                    <td>Samples Collected</td>
                </tr>
                <% ctr = 0 %>
                %for collector in collectors:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${collector[1]}</td>
                        <td>${collector[2]}</td>
                        <td>${collector[3]}</td>
                        <td>${collector[4]}</td>
                        <td><a href="${h.url_for(controller='samples', action='collected_by', sort_id='default', order='default', collector_id=collector[0], last_name=collector[1], first_name=collector[2])}">Samples</a></td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

